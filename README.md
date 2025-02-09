# RNA Structure Evaluation

This project evaluates RNA secondary structures by estimating Gibbs free energy.  It uses a simplified approach based on interatomic distance frequencies and pre-calculated scoring profiles.

## Overview

The project consists of three main parts:

1. **Training:** Generates scoring profiles from known PDB structures.
2. **Evaluation:** Estimates Gibbs free energy for predicted structures.
3. **Plotting:** Visualizes the scoring profiles.

## Usage

1. **Training:**
```bash
   python training.py
```
Configure **RNA_pdb_folder** and **output_folder** in training.py.

Evaluation:
```Bash

python evaluate.py
```
Set pdb_file and scores_dir in evaluate.py.

Plotting:
```Bash

    python plot_profiles.py
```
    Set scores_directory and output_directory in plot_profiles.py.

Files

    training.py: Generates scoring profiles.
    evaluate.py: Evaluates structures.
    plot_profiles.py: Creates plots.
    pdb_processing.py: Reads PDB files.

Data

    PDB files should be placed in the directory specified by RNA_pdb_folder in training.py.
    Scoring profiles (output from training.py) should be used as input for evaluate.py and plot_profiles.py.

Requirements

    Python 3
    matplotlib (pip install matplotlib)

Author

Your Name - Your Email
