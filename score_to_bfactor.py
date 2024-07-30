#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore")


from Bio.PDB import *
import pandas as pd
import sys

'''Script for writing rosetta scores from a scored pdb 
file into the bfactor column of the Calpha atoms. \n
USAGE: score_to_bfactor.py input.pdb output.pdb scoretype \n
Scoretype must be a column in the residue scores at the bottom if the pdb file.'''

if len(sys.argv) < 4:
    print(f"USAGE: {sys.argv[0]} input.pdb output.pdb scoretype")

file = sys.argv[1]
output = sys.argv[2]
scoretype = sys.argv[3]

#file = ""
#output=""
#find out in which line the residue scores start and end
with open(file) as f:
    for i, line in enumerate(f.readlines()):
        if '#BEGIN_POSE_ENERGIES_TABLE' in line: header = i+1
        
        if 'pose' in line: start = i
        if '#END_POSE_ENERGIES_TABLE' in line: end = i -1
rows_toskip = i -end

#read scores as dataframe; skip lines without residue scores
df = pd.read_csv(file, delim_whitespace=True, header=header, skiprows=range(header+1,start+1), skipfooter=rows_toskip)
pd.set_option('display.max_columns', None)
print("these are the residue scores: \n", df, f"\n writing scoretype {scoretype} into bfactor")
#open structure, change bfactor for CA (except for molecules without CA) and save as pdb
p = PDBParser() 
structure = p.get_structure("X", file)[0]
score_list = df[scoretype].tolist()
for score, residue in zip(score_list, structure.get_residues()):
    try:
        atom = residue["CA"]
        atom.set_bfactor(score)
    except:
        for atom in residue:
            atom.set_bfactor(score)

io = PDBIO()
io.set_structure(structure)
io.save(output)

print(f'output saved to {output}. Color in pymol with e.g. the command: spectrum b, blue_white_red, minimum=-10, maximum=10')
#df.to_csv("residue_scores.csv", sep=','); you can also save scores to csv here 

