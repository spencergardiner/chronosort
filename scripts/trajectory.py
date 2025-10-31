# src/trajectory.py
import os
import numpy as np
import math
from Bio.PDB import MMCIFParser, Superimposer

def get_backbone_atoms(model):
    backbone = ['CA']
    return [
        a for a in model.get_atoms()
        if a.get_name() in backbone and a.get_parent().get_id()[0] == ' '
    ]

def compute_rmsd(model1, model2):
    atoms1 = get_backbone_atoms(model1)
    atoms2 = get_backbone_atoms(model2)
    if len(atoms1) != len(atoms2):
        raise ValueError("Structures have different numbers of backbone atoms")
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    return sup.rms

def find_min_rmsd_path(dist):
    n = dist.shape[0]
    best_sum = math.inf
    best_path = None
    for start in range(n):
        path = [start]
        visited = {start}
        current = start
        total = 0.0
        while len(path) < n:
            next_node = min(
                (j for j in range(n) if j not in visited),
                key=lambda j: dist[current, j],
                default=None
            )
            if next_node is None:
                break
            total += dist[current, next_node]
            path.append(next_node)
            visited.add(next_node)
            current = next_node
        if total < best_sum:
            best_sum = total
            best_path = path
    return best_path, best_sum

def make_traj(cif_dir, output_file="ordered.pdb"):
    parser = MMCIFParser(QUIET=True)
    files = [f for f in os.listdir(cif_dir) if f.endswith(".cif")]
    structs = []
    for filename in files:
        filepath = os.path.join(cif_dir, filename)
        try:
            structure = parser.get_structure(filename, filepath)
            model = structure[0]  # Assume first model
            structs.append(model)
        except Exception as e:
            print(f"Error parsing {filename}: {e}")

    n = len(structs)
    if n == 0:
        raise RuntimeError("No valid .cif files found in the specified directory.")
    else:
        dist = np.full((n, n), math.inf)
        np.fill_diagonal(dist, 0.0)
        
        for i in range(n):
            for j in range(i + 1, n):
                try:
                    rmsd = compute_rmsd(structs[i], structs[j])
                    dist[i, j] = dist[j, i] = rmsd
                except ValueError as e:
                    print(f"Error computing RMSD between {files[i]} and {files[j]}: {e}")
        
        best_path, best_sum = find_min_rmsd_path(dist)
        if not best_path: 
            raise RuntimeError("Failed to find an optimal ordering of structures.")
        
        if best_path:
            ordered_files = [files[idx] for idx in best_path]
            print("Optimal ordering (filenames):")
            for filename in ordered_files:
                print(filename)
            print(f"\nTotal sum of RMSDs: {best_sum:.4f}")
     
            # Copy models and superimpose sequentially
            ordered_models = [structs[idx].copy() for idx in best_path]
            for i in range(1, len(ordered_models)):
                prev = ordered_models[i-1]
                curr = ordered_models[i]
                atoms_prev = get_backbone_atoms(prev)
                atoms_curr = get_backbone_atoms(curr)
                sup = Superimposer()
                sup.set_atoms(atoms_prev, atoms_curr)
                sup.apply(curr.get_atoms())  # Apply to all atoms in current model
            
            # Save as multi-model PDB by modifying original files' coordinates
            with open(output_file, 'w') as out:
                for k, idx in enumerate(best_path):
                    filepath = os.path.join(cif_dir, files[idx]) 
                    with open(filepath, 'r') as f:
                        lines = f.readlines()
                    out.write(f"MODEL      {k+1}\n")
                    atom_idx = 0
                    model_atoms = list(ordered_models[k].get_atoms())
                    for line in lines:
                        if not line.startswith('ATOM  '):
                            continue  # skip HETATM and header/footer lines
                        x, y, z = model_atoms[atom_idx].coord
                        atom = model_atoms[atom_idx]
                        resname = atom.get_parent().get_resname()
                        resid = atom.get_parent().get_id()[1]
                        name = atom.get_name()
                        chain_id = atom.get_parent().get_parent().id
                        new_line = f"ATOM  {atom_idx+1:5d} {name:^4} {resname:>3} {chain_id:1}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
                        out.write(new_line)
                        atom_idx += 1
                    out.write("ENDMDL\n")
            print(f"\nSaved ordered structures to '{output_file}' (multi-model PDB)")