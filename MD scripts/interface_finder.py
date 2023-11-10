import sys
import numpy as np
import scipy as sp
import pandas as pd
from Bio.PDB import *
from Bio import BiopythonWarning
import warnings

warnings.simplefilter('ignore', BiopythonWarning)

parser = PDBParser()

pdb_list = open(sys.argv[1]).read().splitlines()
distance_table = pd.DataFrame(columns = ['Res_A', 'Res_B', 'Distance', 'Frame'])

for item in pdb_list:
	pdb_file = open(item)
	structure = parser.get_structure("complex", pdb_file)

	chain_A = structure[0]['A']
	chain_B = structure[0]['B']
	
	for residue_A in chain_A:
		for residue_B in chain_B:
			if residue_A.resname == "GLY":
				atom_A = "CA"
			else:
				atom_A = "CB"
			if residue_B.resname == "GLY":
				atom_B = "CA"
			else:
				atom_B = "CB"	
			distance = residue_A[atom_A] - residue_B[atom_B]
			if distance < 8.0:
				#print(f"A:{residue_A.resname}{residue_A.get_id()[1]},B:{residue_B.resname}{residue_B.get_id()[1]},{atom_B}, Distance:{distance}")
				new_row = {
					'Res_A':f"{residue_A.resname}{residue_A.get_id()[1]}", 
					'Res_B':f"{residue_B.resname}{residue_B.get_id()[1]}",
					'Distance':distance,
					'Frame':pdb_file.name
					}
				new_frame = pd.DataFrame(new_row, index = [0])
				distance_table = pd.concat([distance_table, new_frame], ignore_index = True)
	print(pdb_file.name)

distance_table.to_csv('distance table.csv')
print(distance_table)
