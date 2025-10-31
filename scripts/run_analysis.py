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
    parser.add_argument("--cif_dir", required=True, help="Directory containing .cif files")
    parser.add_argument("--trajectory_file", default="output/trajectory.pdb", help="Output trajectory file")
    parser.add_argument("--vecs_file", default="output/vecs.txt", help="Output eigenvectors file")
    parser.add_argument("--projection_file", default="output/projection.pdb", help="Output PCA projection file")
    parser.add_argument("--scale", type=float, default=30.0, help="Scale for PCA projection")
    parser.add_argument(
        "--components",
        nargs="+",
        type=int,
        default=[0],
        help="One or more 0-based PCA component indices to include in the projection",
    )

    args = parser.parse_args()

    # Convert cif_dir to an absolute path so downstream modules see the same location
    if not os.path.isabs(args.cif_dir):
        args.cif_dir = os.path.abspath(os.getcwd(), args.cif_dir)

    # make sure output directories exist (minimal defensive change)
    _ensure_parent_dir(args.trajectory_file)
    _ensure_parent_dir(args.vecs_file)
    _ensure_parent_dir(args.projection_file)

    # Step 1: build trajectory from CIFs
    make_traj(args.cif_dir, args.trajectory_file)

    # ---- IMPORTANT: call pca_projection with arguments in the order your module expects ----
    # components selects which PCA axes (0-based indices) are combined when generating the projection trajectory
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
