import os
import math
from collections import defaultdict
from training import extract_c3_atoms , compute_distances  # Import from pdb_processing.py

# Constants
MAX_SCORE = 10

def interpolate_score(base1, base2, distance, scores):
    """Find the interpolated score."""
    pair_type = ''.join(sorted([base1, base2]))
    if pair_type not in scores:
        return 0

    score_values = scores[pair_type]
    bin_index = min(int(distance), len(score_values) - 1)  # Closest bin
    score = score_values[bin_index]
    return score

def evaluate_structure(pdb_file, scores):  # scores is now an argument
    """Evaluate an RNA structure."""
    print(f"Evaluating structure: {pdb_file}")

    c3_atoms = extract_c3_atoms(pdb_file)
    if not c3_atoms:
        print(f"Error: No C3' atoms found in {pdb_file}")
        return None

    total_score = 0
    interaction_count = 0

    distances = compute_distances(c3_atoms)  # Use compute_distances here

    for base1, base2, distance in distances:  # Iterate over computed distances
        energy = interpolate_score(base1, base2, distance, scores)
        total_score += energy
        interaction_count += 1

    print(f"Evaluation complete. Total interactions: {interaction_count}")
    print(f"Estimated Gibbs Free Energy: {total_score:.4f}")
    return total_score

def load_scores(scores_dir):
    """Loads scores from files."""
    scores = {}
    for filename in os.listdir(scores_dir):
        if filename.endswith("_scores.txt"):
            base_pair = filename.replace("_scores.txt", "")
            filepath = os.path.join(scores_dir, filename)
            try:
                with open(filepath, 'r') as f:
                    score_values = [float(line.strip()) for line in f]
                    scores[base_pair] = score_values
            except FileNotFoundError:
                print(f"Warning: Score file not found for {base_pair}")
                scores[base_pair] = [MAX_SCORE] * 20  # Default score
    return scores


if __name__ == "__main__":
    pdb_data_dir = "/Users/abdou/Desktop/m2/tp/Main_project/PDB_data/"  # Directory name
    scores_dir = "/Users/abdou/Desktop/m2/tp/Main_project/1st_res"  

    scores = load_scores(scores_dir)
    if not scores:
        print("Error: No scores data found.")
        exit()

    for filename in os.listdir(pdb_data_dir):
        if filename.endswith(".pdb"):  # Filter for .pdb files
            pdb_file = os.path.join(pdb_data_dir, filename)  # Create full path
            evaluate_structure(pdb_file, scores)
