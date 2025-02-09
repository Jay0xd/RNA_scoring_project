import os
import math
from collections import defaultdict

# Constants
BASE_PAIRS = {"AA", "AU", "AC", "AG", "UU", "UC", "UG", "CC", "CG", "GG"}
DISTANCE_BINS = [(i, i + 1) for i in range(20)]
MAX_SCORE = 10
MIN_FREQUENCY = 1e-10

def parse_c3_prime(line):
    if line.startswith("ATOM") and "C3'" in line[12:16]:
        parts = line.split()  # Split the line into parts
        try:
            residue_name = parts[3]  # Example: Residue name
            residue_number = int(parts[5])  # Example: Residue number (integer)
            chain_id = parts[4]  # Example: Chain ID
            x_coord = float(parts[6])  # Example: X coordinate (float)
            y_coord = float(parts[7])  # Example: Y coordinate (float)
            z_coord = float(parts[8])  # Example: Z coordinate (float)
            return residue_name, residue_number, chain_id, x_coord, y_coord, z_coord
        except (ValueError, IndexError):  # Handle potential errors
            return None  # Or raise the exception if you prefer
    return None

def extract_c3_atoms(pdb_file):
    c3_atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            atom_data = parse_c3_prime(line)
            if atom_data:
                c3_atoms.append(atom_data)
    return c3_atoms

def compute_distances(atoms):
    distances = []
    n = len(atoms)
    for i in range(n - 3):
        for j in range(i + 4, n):
            if atoms[i][2] == atoms[j][2]:
                dist = math.sqrt(sum((atoms[i][k] - atoms[j][k])**2 for k in range(3, 6)))
                distances.append((atoms[i][0], atoms[j][0], dist))
    return distances

def calculate_frequencies(distances):
    observed = defaultdict(lambda: [0] * len(DISTANCE_BINS))
    reference = [0] * len(DISTANCE_BINS)

    for res1, res2, dist in distances:
        pair = "".join(sorted((res1, res2)))
        for i, (lower, upper) in enumerate(DISTANCE_BINS):
            if lower <= dist < upper:
                if pair in BASE_PAIRS:
                    observed[pair][i] += 1
                reference[i] += 1
                break

    observed_freq = {}
    for bp, counts in observed.items():
        total_count = sum(counts)
        observed_freq[bp] = [(count / total_count if total_count > 0 else MIN_FREQUENCY) for count in counts]

    total_reference = sum(reference)
    reference_freq = [(count / total_reference if total_reference > 0 else MIN_FREQUENCY) for count in reference]

    return observed_freq, reference_freq

def calculate_scores(observed_freq, reference_freq):
    scores = {}
    for bp, obs_freq_list in observed_freq.items():
        scores[bp] = []
        for obs_freq, ref_freq in zip(obs_freq_list, reference_freq):
            score = -math.log(obs_freq / ref_freq) if ref_freq > 0 and obs_freq > 0 else MAX_SCORE
            scores[bp].append(min(score, MAX_SCORE))
    return scores

def save_scores(scores, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for bp, score_values in scores.items():
        filename = os.path.join(output_dir, f"{bp}_scores.txt")
        with open(filename, 'w') as f:
            for value in score_values:
                f.write(f"{value:.4f}\n")

def main(pdb_folder, output_folder):
    for pdb_file in os.listdir(pdb_folder):
        if pdb_file.endswith(".pdb"):
            pdb_path = os.path.join(pdb_folder, pdb_file)
            c3_atoms = extract_c3_atoms(pdb_path)
            distances = compute_distances(c3_atoms)
            observed_freq, reference_freq = calculate_frequencies(distances)
            scores = calculate_scores(observed_freq, reference_freq)
            save_scores(scores, output_folder)
            print(f"Scores calculated and saved for {pdb_file}")

if __name__ == "__main__":
    rna_pdb_folder = "/Users/abdou/Desktop/m2/tp/Main_project/PDB_data/"
    output_folder = "/Users/abdou/Desktop/m2/tp/Main_project/1st_res/"
    main(rna_pdb_folder, output_folder)

