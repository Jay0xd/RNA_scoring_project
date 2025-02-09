import math
import os
from collections import defaultdict

# Constants
BASE_PAIRS = ["AA", "AU", "AC", "AG", "UU", "CU", "GU", "CC", "CG", "GG"]
DISTANCE_BINS = [(i, i + 1) for i in range(20)]
MAX_SCORE = 10
MIN_FREQUENCY = 1e-10

# Initialize dictionaries
observed_counts = {bp: [0] * len(DISTANCE_BINS) for bp in BASE_PAIRS}
reference_counts = [0] * len(DISTANCE_BINS)

def parse_c3_prime(line):
    if line.startswith("ATOM") and "C3'" in line[12:16]:
        parts = line.split()
        try:
            residue_name = parts[3]
            chain_id = parts[4]
            x_coord = float(parts[6])
            y_coord = float(parts[7])
            z_coord = float(parts[8])
            return residue_name, chain_id, x_coord, y_coord, z_coord
        except (ValueError, IndexError):
            return None
    return None

def extract_c3_atoms(pdb_file):
    c3_atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            atom_data = parse_c3_prime(line)
            if atom_data:
                c3_atoms.append(atom_data)
    return c3_atoms

def compute_distances(c3_atoms):
    n = len(c3_atoms)
    for i in range(n):
        for j in range(i + 4, n):
            if c3_atoms[i][1] == c3_atoms[j][1]:  # Check atoms are in the same chain
                res1, res2 = c3_atoms[i], c3_atoms[j]
                distance = math.sqrt((res1[2] - res2[2])**2 + (res1[3] - res2[3])**2 + (res1[4] - res2[4])**2)
                yield res1[0], res2[0], distance

def update_counts(base1, base2, distance):
    for i, (lower, upper) in enumerate(DISTANCE_BINS):
        if lower <= distance < upper:
            pair_type = ''.join(sorted([base1, base2]))
            if pair_type in BASE_PAIRS:
                observed_counts[pair_type][i] += 1
            reference_counts[i] += 1
            break

def compute_frequencies():
    observed_freqs = {}
    for bp, counts in observed_counts.items():
        total_count = sum(counts)
        observed_freqs[bp] = [(count / total_count if total_count > 0 else MIN_FREQUENCY) for count in counts]

    total_reference = sum(reference_counts)
    reference_freqs = [(count / total_reference if total_reference > 0 else MIN_FREQUENCY) for count in reference_counts]

    return observed_freqs, reference_freqs

def compute_scores(observed_freqs, reference_freqs):
    scores = {}
    for bp, obs_freq_list in observed_freqs.items():
        scores[bp] = []
        for obs_freq, ref_freq in zip(obs_freq_list, reference_freqs):
            safe_ref_freq = ref_freq if ref_freq > 0 else MIN_FREQUENCY
            safe_obs_freq = obs_freq if obs_freq > 0 else MIN_FREQUENCY
            score = -math.log(safe_obs_freq / safe_ref_freq)
            scores[bp].append(min(score, MAX_SCORE))
    return scores

def save_scores(scores, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for bp, score_values in scores.items():
        filename = os.path.join(output_dir, f"{bp}_scores")
        with open(filename, 'w') as f:
            for value in score_values:
                f.write(f"{value:.4f}\n")

def main(pdb_folder, output_folder):
    for pdb_file in os.listdir(pdb_folder):
        if pdb_file.endswith(".pdb"):
            pdb_path = os.path.join(pdb_folder, pdb_file)
            c3_atoms = extract_c3_atoms(pdb_path)

            for base1, base2, distance in compute_distances(c3_atoms):
                update_counts(base1, base2, distance)

            observed_freqs, reference_freqs = compute_frequencies()
            scores = compute_scores(observed_freqs, reference_freqs)
            save_scores(scores, output_folder)
            print(f"Scores calculated and saved for {pdb_file}")


if __name__ == "__main__":
    rna_pdb_folder = "/Users/abdou/Desktop/m2/tp/Main_project/PDB_data/"  # Corrected path
    output_folder = "/Users/abdou/Desktop/m2/tp/Main_project/score_output_profiles/"  # Corrected path
    main(rna_pdb_folder, output_folder)