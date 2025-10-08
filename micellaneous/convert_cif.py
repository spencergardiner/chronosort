#!/usr/bin/env python3
import sys
import os
from Bio.PDB import MMCIFParser, PDBIO

def convert_cif_to_pdb(cif_file, pdb_file):
    parser = MMCIFParser()
    try:
        structure = parser.get_structure('structure', cif_file)
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_file)
        return True
    except Exception as e:
        print(f"Error converting {cif_file}: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_cif.py <cif_file> <pdb_file>")
        sys.exit(1)
    
    cif_file = sys.argv[1]
    pdb_file = sys.argv[2]
    
    success = convert_cif_to_pdb(cif_file, pdb_file)
    sys.exit(0 if success else 1)
