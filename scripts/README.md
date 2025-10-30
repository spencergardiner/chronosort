# Protein Ensemble PCA Analysis

This project builds a trajectory from a set of `.cif` protein structures, computes PCA on the aligned backbone atoms, and generates a projected trajectory along principal components.

## Features

- RMSD-based optimal ordering of structures
- Multi-model PDB trajectory generation
- PCA projection with customizable components and scale
- Outputs eigenvectors, projected trajectory, and eigenvalue plot

## Getting Started

First, clone the repository and navigate into it: 

```bash
git clone https://github.com/dellacortelab/chronosort.git
cd chronosort
```

## Requirements

- Python 3.8+
- [MDAnalysis](https://www.mdanalysis.org/)
- [Biopython](https://biopython.org/)
- numpy
- matplotlib

Next, install dependencies with:

``` bash
pip install -r scripts/requirements.txt
```

## Usage

If you'd like to verify your setup using the provided test sequences, run:

```bash
python scripts/run_analysis.py --cif_dir test_data/
```

Confirm that `test_data/` contains `.cif` files as expected. This will generate all output files in the `output/` directory.


Once your own `.cif` files are ready, run the pipeline with:

```bash
python scripts/run_analysis.py --cif_dir path/to/cifs
```

### What it does:
- Constructs a trajectory from CIF files using RMSD-based optimal ordering
- Aligns backbone atoms and performs PCA
- Outputs:
  - `trajectory.pdb`: ordered multi-model trajectory
  - `projection.pdb`: PCA-projected trajectory
  - `vecs.txt`: saved eigenvectors
  - `eigenvalues.png`: scree plot of PCA variance

### Modules involved:
- `trajectory.py`: builds the trajectory and superimposes structures
- `pca_analysis.py`: performs PCA and generates projections

### Optional arguments:
Customize the analysis with these flags:

```bash
--trajectory_file     Output trajectory file (default: output/trajectory.pdb)
--vecs_file           Output eigenvectors file (default: output/vecs.txt)
--projection_file     Output PCA projection file (default: output/projection.pdb)
--scale               Scale for PCA projection (default: 30.0)
--components          Number of PCA components to use (default: [0])
```