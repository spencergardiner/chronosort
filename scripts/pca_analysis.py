# src/pca_analysis.py
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align, pca
# from source.trajectory import make_traj  # No longer needed

def pca_projection(path_to_pdbs, trajectory_file, vecs_file, projection_file, sel_str="name CA", components=[0], projection=[0, 1], scale=30.0):
    evec=[]
    # make_traj(path_to_pdbs, trajectory_file)  # Removed to prevent duplicate printing
    u = mda.Universe(trajectory_file)
    sel = u.select_atoms(sel_str)
    print(f"Selected {sel.n_atoms} CA atoms from trajectory.")
    if sel.n_atoms == 0:
        raise ValueError(f"No atoms selected with selection string '{sel_str}'. Check your input structure.")
    align.AlignTraj(u, u, select=sel_str, in_memory=True).run()
    pc = pca.PCA(u, select=sel_str, align=True)
    pc.run()

    n_atoms = sel.atoms.n_atoms

    frame0 = sel.atoms.positions.copy()
    for n in projection:
        evec.append(pc.results.p_components[:, n].reshape(n_atoms, 3))
    scales = np.linspace(-scale, scale, 201)
    combined_vector=np.zeros(n_atoms * 3).reshape(n_atoms, 3)
    for i in evec:
        combined_vector += i
        with open(vecs_file, "a") as f: # this might need to be changed to a different file
            np.savetxt(f, i)
    with open(vecs_file, "w") as f:   # open once in write mode
        for idx, n in enumerate(projection, start=1):
            vec = pc.results.p_components[:, n].reshape(n_atoms, 3)
            evec.append(vec)
            eigval = pc.results.variance[n]   # correct for older MDAnalysis

            f.write(f"PC{idx} eigenvalue={eigval:.6f} total={np.sum(pc.results.variance)}\n")
            np.savetxt(f, vec)
            f.write("\n")

    variance = pc.results.variance
    percent_variance = 100 * variance / np.sum(variance)
    plt.figure(figsize=(8, 5))
    plt.plot(percent_variance[:30], marker='o')
    plt.xlabel("Principal Component Index")
    plt.ylabel("Percent of Variance Explained")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("output/eigenvalues.png")  # optional: save to output folder
    plt.close()
    print(f"Saved eigenvalues to 'output/eigenvalues.png'")


    frames = [frame0 + s * combined_vector for s in scales]
    with mda.Writer(projection_file, n_atoms=n_atoms) as W:
        for coords in frames:
            sel.atoms.positions = coords
            W.write(sel.atoms)
    print(f"Saved PCA projection to '{projection_file}'")