# scripts/run_analysis.py
import argparse
import os
import warnings
from trajectory import make_traj
from pca_analysis import pca_projection

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

def _ensure_parent_dir(path):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)

def main():
    parser = argparse.ArgumentParser(description="Protein PCA Analysis")
    parser.add_argument("--cif_dir", required=True, help="Directory with .cif files")
    parser.add_argument("--trajectory_file", default="output/trajectory.pdb", help="Output trajectory file")
    parser.add_argument("--vecs_file", default="output/vecs.txt", help="Output eigenvectors file")
    parser.add_argument("--projection_file", default="output/projection.pdb", help="Output PCA projection file")
    parser.add_argument("--scale", type=float, default=30.0, help="Scale for PCA projection")
    parser.add_argument("--components", nargs="+", type=int, default=[0], help="PCA components to use")

    args = parser.parse_args()

    import os

    # Convert cif_dir to an absolute path relative to the script's location if it's relative
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    args.cif_dir = (
        args.cif_dir
        if os.path.isabs(args.cif_dir)
        else os.path.join(SCRIPT_DIR, args.cif_dir)
    )

    # make sure output directories exist (minimal defensive change)
    _ensure_parent_dir(args.trajectory_file)
    _ensure_parent_dir(args.vecs_file)
    _ensure_parent_dir(args.projection_file)

    # Step 1: build trajectory from CIFs
    make_traj(args.cif_dir, args.trajectory_file)

    # ---- IMPORTANT: call pca_projection with arguments in the order your module expects ----
    # pca_projection(path_to_pdbs, trajectory_file, vecs_file, projection_file, sel_str=..., components=..., projection=..., scale=...)
    pca_projection(
        args.cif_dir,
        args.trajectory_file,
        args.vecs_file,
        args.projection_file,
        components=args.components,
        scale=args.scale
    )

if __name__ == "__main__":
    main()
