#!/bin/bash

# This script:
# 1. Concatenates all CIF files in a folder into a GROMACS trajectory
# 2. Calculates pairwise CA RMSD between all structures
# 3. Creates a new trajectory with frames sorted by RMSD similarity

# Usage: ./sort_cif_by_rmsd.sh <cif_folder> <output_prefix>
rm *# 
# Check if proper arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <cif_folder> <output_prefix>"
    exit 1
fi

# Assign arguments to variables
CIF_FOLDER="$1"
OUTPUT_PREFIX="$2"

# Count CIF files
CIF_COUNT=$(ls -1 "$CIF_FOLDER"/*.cif 2>/dev/null | wc -l)
if [ "$CIF_COUNT" -eq 0 ]; then
    echo "Error: No CIF files found in $CIF_FOLDER"
    exit 1
fi

echo "Found $CIF_COUNT CIF files in $CIF_FOLDER"

# Check for required dependencies
command -v gmx >/dev/null 2>&1 || { echo "Error: GROMACS is required but not installed. Exiting."; exit 1; }
python3 -c "import Bio" >/dev/null 2>&1 || { echo "Error: Biopython is required but not installed. Install with 'pip install biopython'. Exiting."; exit 1; }

# Create a temporary directory for intermediate files
TEMP_DIR=$(mktemp -d)
echo "Created temporary directory: $TEMP_DIR"

# Step 1: Prepare and concatenate all CIF files into a single trajectory
echo "Step 1: Preparing and concatenating CIF files into a trajectory..."

# Select first CIF as reference for topology
FIRST_CIF=$(ls -1 "$CIF_FOLDER"/*.cif | head -n 1)
echo "Using $FIRST_CIF as reference for topology"

# First, convert CIF to PDB (GROMACS pdb2gmx works better with PDB format)
echo "Converting reference CIF to PDB format..."
mkdir -p "$TEMP_DIR/converted"
python3 -c "
import sys
try:
    from Bio.PDB import MMCIFParser, PDBIO
    parser = MMCIFParser()
    structure = parser.get_structure('reference', '$FIRST_CIF')
    io = PDBIO()
    io.set_structure(structure)
    io.save('$TEMP_DIR/converted/reference.pdb')
    print('Conversion successful')
except ImportError:
    print('Error: Biopython is required for CIF conversion')
    sys.exit(1)
except Exception as e:
    print(f'Error during conversion: {e}')
    sys.exit(1)
"

# Check if conversion was successful
if [ ! -f "$TEMP_DIR/converted/reference.pdb" ]; then
    echo "Error: Failed to convert CIF to PDB. Ensure Biopython is installed (pip install biopython)"
    exit 1
fi

# Create a topology file from the converted PDB
gmx pdb2gmx -f "$TEMP_DIR/converted/reference.pdb" -o "$TEMP_DIR/reference.gro" -p "$TEMP_DIR/topol.top" -water none -ff amber99sb-ildn -ignh <<EOF
1
EOF

# Create a simulation box
gmx editconf -f "$TEMP_DIR/reference.gro" -o "$TEMP_DIR/reference_box.gro" -c -d 1.0 -bt cubic


# Create an index file to select CA atoms
echo "[ CA ]" > "$TEMP_DIR/ca_index.ndx"
echo "name CA" >> "$TEMP_DIR/ca_index.ndx"

# Process each CIF file and create temporary trajectory parts
mkdir -p "$TEMP_DIR/parts"
mkdir -p "$TEMP_DIR/converted"
COUNT=0


for CIF_FILE in "$CIF_FOLDER"/*.cif; do
    BASENAME=$(basename "$CIF_FILE")
    echo "Processing $BASENAME"
    
    # Convert CIF to PDB
    CONVERTED_PDB="$TEMP_DIR/converted/$COUNT.pdb"
    python3 "./convert_cif.py" "$CIF_FILE" "$CONVERTED_PDB"
    
    if [ ! -f "$CONVERTED_PDB" ]; then
        echo "Warning: Failed to convert $BASENAME. Skipping."
        continue
    fi
    
    # Process the converted PDB with GROMACS
    # (This assumes all structures have the same protein topology)
    gmx editconf -f "$CONVERTED_PDB" -o "$TEMP_DIR/parts/$COUNT.gro" -resnr 1
    
    COUNT=$((COUNT + 1))
done

#Here, I need to have the PDB files sorted into the right order

cat "$TEMP_DIR/converted/*pdb" > "$TEMP_DIR/traj.pdb"

gmx pdb2gmx -f "$TEMP_DIR/converted/0.pdb" -o "$TEMP_DIR/conf.gro" -p "$TEMP_DIR/topol.top"<<EOF
6
1
EOF

gmx editconf -f "$TEMP_DIR/conf.gro" -o "$TEMP_DIR/boxed.gro" -c -d 1.0 -bt cubic


gmx solvate -cp "$TEMP_DIR/boxed.gro" -cs -o "$TEMP_DIR/solvated.gro" -p "$TEMP_DIR/topol.top"

gmx grompp -f "./ions.mdp" -c "$TEMP_DIR/solvated.gro" -p "$TEMP_DIR/topol.top" -o "$TEMP_DIR/dummy.tpr" -maxwarn 100

python pdb2gmx.py "$TEMP_DIR"

echo $TEMP_DIR
echo "Script completed successfully."