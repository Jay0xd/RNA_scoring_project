import os
import matplotlib.pyplot as plt

def read_scores(directory):
    """Reads score files and returns data as a dictionary."""
    all_scores = {}
    for filename in os.listdir(directory):
        if filename.endswith("_scores.txt"):
            base_pair = filename.replace("_scores.txt", "")
            filepath = os.path.join(directory, filename)
            try:
                with open(filepath, 'r') as f:
                    scores = [float(line.strip()) for line in f]
                    all_scores[base_pair] = scores
            except FileNotFoundError:
                print(f"Warning: Score file not found for {base_pair}")
    return all_scores

