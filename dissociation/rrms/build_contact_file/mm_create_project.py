#!python
import os
import argparse
import sys
from time import sleep
import subprocess
import fileinput
import platform
import open3SPN2

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.insert(0, OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

from openmmawsem import *
from helperFunctions.myFunctions import *
import ffAWSEM

parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create \
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("protein", help="The name of the protein, \
            do: python3 ~/OPENAWSEM_LOCATION/mm_create_project.py")

args = parser.parse_args()
proteinName = pdb_id = args.protein
pdb = f"{pdb_id}.pdb"



#Fix the system (adds missing atoms)
fix=open3SPN2.fixPDB(pdb)
#Create a table containing both the proteins and the DNA
complex_table=open3SPN2.pdb2table(fix)
#Coarse Grain the system
dna_atoms=open3SPN2.DNA.CoarseGrain(complex_table)
protein_atoms=ffAWSEM.Protein.CoarseGrain(complex_table)
#Merge the models
import pandas
Coarse=pandas.concat([dna_atoms,protein_atoms],sort=False)
Coarse.index=range(len(Coarse))
Coarse.serial=list(Coarse.index)
#Save the protein_sequence
from Bio.PDB.Polypeptide import three_to_one
_AWSEMresidues=['IPR','IGL','NGP']
protein_data=Coarse[Coarse.resname.isin(_AWSEMresidues)].copy()
resix = (protein_data.chainID + '_' + protein_data.resSeq.astype(str))
res_unique = resix.unique()
protein_data['resID'] = resix.replace(dict(zip(res_unique, range(len(res_unique)))))
protein_sequence=[r.iloc[0]['real_resname'] for i, r in protein_data.groupby('resID')]
protein_sequence_one = [three_to_one(a) for a in protein_sequence]

with open('protein.seq','w+') as ps:
    ps.write(''.join(protein_sequence_one))
# Create a merged PDB
def writePDB(atoms,pdb_file):
    with open(pdb_file, 'w+') as pdb:
        for i, atom in atoms.iterrows():
            pdb_line = f'{atom.recname:<6}{atom.serial:>5} {atom["name"]:^4}{atom.altLoc:1}'+\
                       f'{atom.resname:<3} {atom.chainID:1}{atom.resSeq:>4}{atom.iCode:1}   '+\
                       f'{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}' +\
                       f'{atom.occupancy:>6.2f}{atom.occupancy:>6.2f}'+' ' * 10 +\
                       f'{atom.element:>2}{atom.charge:>2}'
            assert len(pdb_line) == 80, f'An item in the atom table is longer than expected ({len(pdb_line)})\n{pdb_line}'
            pdb.write(pdb_line + '\n')
writePDB(Coarse,'clean.pdb')



