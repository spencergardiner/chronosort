import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser, Selection, Superimposer
from pathlib import Path
import pandas as pd

def find_closest_neighbors(rmsd_matrix):
    """
    Find the closest neighbor for each index in the distance matrix.
    
    Parameters:
    -----------
    rmsd_matrix : numpy.ndarray
        A square distance matrix
        
    Returns:
    --------
    sorted_list : list
        A list where each element is the index of the closest neighbor
    """
    # Create a copy to avoid modifying the original matrix
    working_matrix = rmsd_matrix.copy()
    
    # Set diagonal elements to NaN (distance to self)
    mask = np.eye(working_matrix.shape[0], dtype=bool)
    working_matrix[mask] = np.nan
    
    sorted_list = []
    
    for idx in range(working_matrix.shape[0]):
        row = working_matrix[idx]
        argmin_value = np.nanargmin(row)  # Use nanargmin to ignore NaN values
        sorted_list.append(argmin_value)
        
        # Mask the found value in all rows to avoid selecting it again
        # Set the entire column of argmin_value to NaN
        working_matrix[:, argmin_value] = np.nan
        
    return sorted_list

def extract_ca_atoms(structure):
    """Extract Carbon Alpha atoms from a structure"""
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    ca_atoms.append(residue["CA"])
    return ca_atoms

def calculate_rmsd(structure1, structure2):
    """Calculate RMSD between two structures using only CA atoms"""
    # Extract CA atoms
    ca_atoms_1 = extract_ca_atoms(structure1)
    ca_atoms_2 = extract_ca_atoms(structure2)
    
    # Check if both structures have the same number of CA atoms
    if len(ca_atoms_1) != len(ca_atoms_2):
        print(f"Warning: Different number of CA atoms: {len(ca_atoms_1)} vs {len(ca_atoms_2)}")
        # Find the minimum length to use for comparison
        min_length = min(len(ca_atoms_1), len(ca_atoms_2))
        ca_atoms_1 = ca_atoms_1[:min_length]
        ca_atoms_2 = ca_atoms_2[:min_length]
    
    # Use Superimposer to calculate RMSD
    super_imposer = Superimposer()
    super_imposer.set_atoms(ca_atoms_1, ca_atoms_2)
    return super_imposer.rms

def main():
    # Check if directory path is provided
    if len(sys.argv) != 2:
        print("Usage: python pairwise_ca_rmsd.py <directory_path>")
        sys.exit(1)
    
    # Get directory path
    dir_path = sys.argv[1]
    
    # Find all PDB files in the directory
    pdb_files = glob.glob(os.path.join(dir_path, "*.pdb"))
    
    if not pdb_files:
        print(f"No PDB files found in {dir_path}")
        sys.exit(1)
    
    print(f"Found {len(pdb_files)} PDB files")
    
    # Initialize PDB parser
    parser = PDBParser(QUIET=True)
    
    # Create a dictionary to store structures
    structures = {}
    for pdb_file in pdb_files:
        file_name = os.path.basename(pdb_file)
        structure_id = os.path.splitext(file_name)[0]
        structures[structure_id] = parser.get_structure(structure_id, pdb_file)
    
    # Calculate pairwise RMSD
    structure_ids = list(structures.keys())
    rmsd_matrix = np.zeros((len(structure_ids), len(structure_ids)))
    
    for i, id1 in enumerate(structure_ids):
        for j, id2 in enumerate(structure_ids):
            if i <= j:  # Only calculate upper triangle (including diagonal)
                rmsd = calculate_rmsd(structures[id1], structures[id2])
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd  # Matrix is symmetric
    
    # Create DataFrame for better visualization
    rmsd_df = pd.DataFrame(rmsd_matrix, index=structure_ids, columns=structure_ids)
    
    # Save results to CSV
    output_csv = os.path.join(dir_path, "pairwise_ca_rmsd.csv")
    rmsd_df.to_csv(output_csv)
    print(f"Saved RMSD matrix to {output_csv}")
    
    # Create heatmap visualization
    plt.figure(figsize=(10, 8))
    sns.heatmap(rmsd_df, annot=True, cmap="viridis", fmt=".3f")
    plt.title("Pairwise CA RMSD between PDB structures")
    plt.tight_layout()
    
    # Save heatmap
    heatmap_path = os.path.join(dir_path, "pairwise_ca_rmsd_heatmap.png")
    plt.savefig(heatmap_path, dpi=300)
    print(f"Saved heatmap to {heatmap_path}")
    
    # Show the heatmap (comment out if running in headless environment)
    # plt.show()


import MDAnalysis as mda
import glob
import argparse

# Step 1: Parse input arguments
parser = argparse.ArgumentParser(description="Process PDB files from a given path.")
parser.add_argument("path", type=str, help="Path to the directory containing PDB files.")
args = parser.parse_args()

# Step 2: Load all PDB files in sorted order
pdb_files = sorted(glob.glob(f"{args.path}/converted/*.pdb"))

# Initialize PDB parser
parser = PDBParser(QUIET=True)

# Create a dictionary to store structures
structures = {}
for pdb_file in pdb_files:
    file_name = os.path.basename(pdb_file)
    structure_id = os.path.splitext(file_name)[0]
    structures[structure_id] = parser.get_structure(structure_id, pdb_file)

# Calculate pairwise RMSD
structure_ids = list(structures.keys())
rmsd_matrix = np.zeros((len(structure_ids), len(structure_ids)))

list = find_closest_neighbors(rmsd_matrix)

# Structure_ids and sorted list just need to be used to resorte the PDB Files
sorted_pdb_files = []
for i in range(len(list)):
    sorted_pdb_files.append(pdb_files[list[i]])

pdf_files = sorted_pdb_files

# Parameters
timestep_ps = 10.0


# Use first structure as reference
u = mda.Universe(pdb_files[0])
with mda.Writer(f"{args.path}/combined.xtc", n_atoms=u.atoms.n_atoms) as W:
    for i, pdb in enumerate(pdb_files):
        temp_u = mda.Universe(pdb)
        ts = temp_u.trajectory[0]
        ts.time = i * timestep_ps
        ts.frame = i
        temp_u.trajectory.ts = ts  # assign modified timestep back
        W.write(temp_u.atoms)
