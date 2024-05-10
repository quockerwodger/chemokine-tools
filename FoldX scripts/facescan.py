# -*- coding: utf-8 -*-

from Bio.PDB import *
from Bio import BiopythonWarning
import warnings

warnings.simplefilter('ignore', BiopythonWarning)

parser = PDBParser()

CKBP_interface = []

#Loading in the pdb
pdb_file = input("Provide the clash pdb: ")
structure = parser.get_structure("complex", pdb_file)
chain_A = structure[0]['A']
chain_C = structure[0]['C']
Chem_interface=[]

#Calculating interface residues
for residue_A in chain_A:
    for residue_B in chain_C:
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
            if 157 <= residue_A.get_id()[1] <= 191:
                pass
                #This exists to exclude a highly flexible loop found within CKBP from calculations
            else:
                CKBP_interface.extend([residue_A])
                Chem_interface.extend([residue_B])


#Removing duplicates from interface list & putting into numerical order
CKBP_interface = set(CKBP_interface)
sorted_CKBP = sorted(CKBP_interface, key=lambda x: int(x.get_id()[1]))
Chem_interface=set(Chem_interface)
print(sorted_CKBP)

# Creating dictionary to convert from 3 letter code to 1 letter code
d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
 'GLY': 'G', 'HIS': 'H', 'HSD':'H', 'HSE':'H', 'HSP':'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

#A dictionary that assigns a number to each amino acid type
aa_dict={1:'A', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I', 9:'K', 10:'L',
         11:'M', 12:'N', 13:'P', 14:'Q', 15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y'}

#Creating the mutant sequences. Will build out sequence based on wt structure until reaches residue from interface list, will then create mutants before starting to build out again
mutant_list=[]
wt_seq=[]
for residue in chain_A:
    if_occured=False
    if residue in CKBP_interface:
        iteration=1
        mutant_seq=[]
        while iteration < 21:
            mutant_seq.extend(wt_seq)
            mutant_seq.extend(aa_dict[iteration])
            #print("poop:\n",''.join(mutant_seq))
            mutant_string=''.join(mutant_seq)
            mutant_list.append(mutant_string)
            iteration += 1
            mutant_seq=[]
        if iteration == 21:
            wt_seq.append(d3to1[residue.resname])
    else:
        mutant_seq=[]
        wt_seq.append(d3to1[residue.resname])

#Combining mutant sequences into single strings & saving output file
wt=[''.join(wt_seq)]
wt.extend(mutant_list)
with open('mutant_file.txt', 'w') as f:
    for line in wt:
        f.write("%s\n" % line)
        

        