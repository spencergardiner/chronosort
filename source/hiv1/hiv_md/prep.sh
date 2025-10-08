#!/bin/bash
#SBATCH --job-name=hemo
#SBATCH --time=72:00:00       # Max runtime
#SBATCH --gres=gpu:1          
#SBATCH --cpus-per-gpu=4
#SBATCH --mem=60G             
#SBATCH --requeue
# module load cuda/12.4.1 gromacs/2024.2-w2p4z73
module load cuda/12.4.1 gromacs/2024.3-qklo76t
export OMP_NUM_THREADS=4

COMPLEX_PDB="$1"

# Check that input file exists
if [[ ! -f "$COMPLEX_PDB" ]]; then
    echo "Error: File '$COMPLEX_PDB' not found."
    exit 1
fi


cp "$COMPLEX_PDB" protein.pdb

printf "1\n1" | gmx pdb2gmx -f protein.pdb -o protein.gro -p topol.top -ignh

gmx editconf -f protein.gro -o newbox.gro -bt dodecahedron -d 1.0 -c

gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

gmx grompp -f mdps/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
echo 13 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
gmx grompp -f mdps/em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -nb gpu -ntmpi 1 -ntomp 4


gmx make_ndx -f em.gro -o index.ndx << 'EOF'
1 | 13
q
EOF

gmx grompp -f mdps/nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt
gmx grompp -f mdps/npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt
gmx grompp -f mdps/md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
gmx mdrun -deffnm md