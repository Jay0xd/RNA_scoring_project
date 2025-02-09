import os
from scoring import read_scores
import matplotlib.pyplot as plt

def plot_interaction_profiles(all_scores, output_dir):
    """Plots interaction profiles."""
    os.makedirs(output_dir, exist_ok=True)

    # General Plot
    print("Generating general interaction profile plot...")
    plt.figure(figsize=(10, 6))
    for base_pair, scores in all_scores.items():
        plt.plot(range(len(scores)), scores, label=base_pair, linewidth=2)  # Use range for x-axis

    plt.title("Interaction Profiles: Score vs Distance", fontsize=16)
    plt.xlabel("Distance (Å)", fontsize=14)
    plt.ylabel("Score", fontsize=14)
    plt.legend(title="Base Pair", fontsize=12)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    general_plot_path = os.path.join(output_dir, "general_interaction_profile.png")
    plt.savefig(general_plot_path, dpi=300)
    plt.close()
    print(f"General plot saved to {general_plot_path}")

    # Individual Plots
    print("Generating individual plots for each base pair...")
    for base_pair, scores in all_scores.items():
        plt.figure(figsize=(8, 5))
        plt.plot(range(len(scores)), scores, label=base_pair, color="b", linewidth=2)

        plt.title(f"Interaction Profile: {base_pair}", fontsize=16)
        plt.xlabel("Distance (Å)", fontsize=14)
        plt.ylabel("Score", fontsize=14)
        plt.legend(fontsize=12)
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.tight_layout()

        individual_plot_path = os.path.join(output_dir, f"{base_pair}_interaction_profile.png")
        plt.savefig(individual_plot_path, dpi=300)
        plt.close()
        print(f"Plot for {base_pair} saved to {individual_plot_path}")

# Example usage:
if __name__ == "__main__":
    scores_directory = "/Users/abdou/Desktop/m2/tp/Main_project/1st_res/"  # Replace with the actual path
    output_directory = "/Users/abdou/Desktop/m2/tp/Main_project/plot_"  # Replace with the desired output path
    all_scores_data = read_scores(scores_directory)
    plot_interaction_profiles(all_scores_data, output_directory)