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
