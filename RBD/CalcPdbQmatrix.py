#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan, David Winogradoff and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys
from VectorAlgebra import *

#from Bio.PDB.PDBParser import PDBParser

atom_type = {'1' : 'C', '2' : 'N', '3' : 'O', '4' : 'C', '5' : 'H', '6' : 'C'}
atom_desc = {'1' : 'C-Alpha', '2' : 'N', '3' : 'O', '4' : 'C-Beta', '5' : 'H-Beta', '6' : 'C-Prime'}
PDB_type = {'1' : 'CA', '2' : 'N', '3' : 'O', '4' : 'CB', '5' : 'HB', '6' : 'C' }

class PDB_Atom:
	no = 0
	ty = ''
	mol = 0
	res = 'UNK'
	res_no = 0
	x = 0.0
	y = 0.0
	z = 0.0
	atm = 'C'
	
	def __init__(self, no, ty, mol, res, res_no, x, y, z, atm):
		self.no = no
		self.ty = ty
		self.mol = mol
		self.res = res
		self.res_no = res_no
		self.x = x
		self.y = y
		self.z = z
		self.atm = atm
		
	def write_(self, f):
		f.write('ATOM')
		f.write(('       '+str(self.no))[-7:])
		f.write('  ')
		f.write(self.mol)
		f.write('  ')
		f.write((self.ty+'    ')[:4])
		f.write(self.res)
		f.write(' ')
		f.write('T')
		f.write(('    '+str(self.res_no))[-4:])
		f.write(('            '+str(round(self.x,3)))[-12:])
		f.write(('        '+str(round(self.y,3)))[-8:])
		f.write(('        '+str(round(self.z,3)))[-8:])
		f.write('  1.00')
		f.write('  0.00')
		f.write(('            '+self.atm)[-12:]+'  ')
		f.write('\n')

class Atom:
	No = 0
	ty = ''
	x = 0.0
	y = 0.0
	z = 0.0
	desc = ''
	
	def __init__(self, No, ty, No_m, x, y, z, desc=''):
		self.No = No
		self.ty = ty
		self.No_m = No_m
		self.x = x
		self.y = y
		self.z = z
		self.desc = desc
	
	def write_(self, f):
		f.write(str(self.No))
		f.write(' ')
		f.write(PDB_type[self.No_m])
		f.write(' ')
		f.write(str(round(self.x,8)))
		f.write(' ')
		f.write(str(round(self.y,8)))
		f.write(' ')
		f.write(str(round(self.z,8)))
		f.write(' ')
		f.write(self.desc)
		f.write('\n')

if len(sys.argv)!=3:                               # and len(sys.argv)!=5:
	print "\nCalcPdbQ.py Output_file qonuchic_flag(1 for q_o, 0 for q_w)\n"
#	print
#	print "\t\t-i\tcalculate individual q values for each chain"
#	print
#	exit()

splitq = False
#for iarg in range(0, len(sys.argv)):
#	if sys.argv[iarg]=="-i":
#		splitq = True
#		sys.argv.pop(iarg)
#
#struct_id = sys.argv[1]
#if struct_id[-4:].lower()==".pdb":
#	pdb_file = struct_id
#else:
#	pdb_file = struct_id + ".pdb"
#
#struct_id2 = sys.argv[2]
#if struct_id2[-4:].lower()==".pdb":
#        pdb_file2 = struct_id2
#else:
#        pdb_file2 = struct_id2 + ".pdb"

output_file = ""
output_file = sys.argv[1]

sigma_exp = 0.15
qo_flag = float(sys.argv[2])

out = open(output_file, 'w')

from Bio.PDB.PDBParser import PDBParser

p = PDBParser(PERMISSIVE=1)

def computeQ():
	if len(ca_atoms_pdb2)!=len(ca_atoms_pdb):
		print "Error. Length mismatch!"
		print "Pdb1: ", len(ca_atoms_pdb), "Pdb2: ", len(ca_atoms_pdb2)
		exit()
	Q = 0.0
	norm = 0
	N = len(ca_atoms_pdb)
	min_sep = 3
        for ia in range(0, 100):
                for ja in range(max(100, ia+min_sep), N):
#	for ia in range(0, N):
#		for ja in range(ia+min_sep, N):
				r = vabs(vector(ca_atoms_pdb[ia], ca_atoms_pdb[ja]))
				rn = vabs(vector(ca_atoms_pdb2[ia], ca_atoms_pdb2[ja]))
                                if qo_flag == 1 and rn >= 9.5 and r >=9.5: continue
                                dr = r - rn
				Q = Q + exp(-dr*dr/(2*sigma_sq[ja-ia]))
				norm = norm + 1
        for ia in range(100, 200):
                for ja in range(max(200, ia+min_sep), N):
                                r = vabs(vector(ca_atoms_pdb[ia], ca_atoms_pdb[ja]))
                                rn = vabs(vector(ca_atoms_pdb2[ia], ca_atoms_pdb2[ja]))
                                if qo_flag == 1 and rn >= 9.5 and r >=9.5: continue
                                dr = r - rn
                                Q = Q + exp(-dr*dr/(2*sigma_sq[ja-ia]))
                                norm = norm + 1
	Q = Q/norm
	return Q

Nframe = 60


for i in range(1,Nframe+1):
     sigma = []
     sigma_sq = []
     ca_atoms_pdb = []
     struct_id = 'end-' + str(i)
     pdb_file = 'end-' + str(i) + '.pdb'
     s = p.get_structure(struct_id, pdb_file)
     chains = s[0].get_list()
     for chain in chains:
     	for res in chain:
     		is_regular_res = res.has_id('CA') and res.has_id('O')
     		res_id = res.get_id()[0]
     	        if is_regular_res:#(res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS' ) and is_regular_res:
     			ca_atoms_pdb.append(res['CA'].get_coord())
     for l in range(0, len(ca_atoms_pdb)+1):
                sigma.append( (1+l)**sigma_exp )
                sigma_sq.append(sigma[-1]*sigma[-1])

     for j in range(1,Nframe+1):
         ca_atoms_pdb2 = []
         struct_id2 = 'end-' + str(j)
         pdb_file2 = 'end-' + str(j) + '.pdb'     
         s2 = p.get_structure(struct_id2, pdb_file2)
         chains = s2[0].get_list()
         for chain in chains:
                 for res in chain:
                         is_regular_res = res.has_id('CA') and res.has_id('O')
                         res_id = res.get_id()[0]
                         if is_regular_res:#(res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS' ) and is_regular_res:
                                 ca_atoms_pdb2.append(res['CA'].get_coord())
         
         if len(ca_atoms_pdb) != len(ca_atoms_pdb2):
         	print "Error: Pdb structures have different lengths!"
         	exit() 
         
         if len(ca_atoms_pdb)>0:
                q = computeQ()
         	out.write(str(round(q,3)))
         	out.write('\t')
     out.write('\n')
out.close()
