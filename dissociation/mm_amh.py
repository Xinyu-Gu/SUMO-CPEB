#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import fileinput
import platform
import importlib.util
import open3SPN2
import numpy as np

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.append(OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

# __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
# __author__ = 'Wei Lu'

from openmmawsem import *
from helperFunctions.myFunctions import *
import ffAWSEM

parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create \
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-c", "--chain", type=str, default="-1")
parser.add_argument("--thread", type=int, default=2, help="default is using 2 CPUs, -1 is using all")
parser.add_argument("-p", "--platform", type=str, default="OpenCL", help="Could be OpenCL, CUDA and CPU")
parser.add_argument("-t", "--trajectory", type=str, default="./movie.pdb")
parser.add_argument("--params", type=str, default="./params.py")
parser.add_argument("-o", "--output", type=str, default=None, help="The location of file that show your energy and Q infomation.")
args = parser.parse_args()

if (args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir


simulation_platform = args.platform
platform = Platform.getPlatformByName(simulation_platform)
if simulation_platform == "CPU":
    if args.thread != -1:
        platform.setPropertyDefaultValue("Threads", str(args.thread))
    print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")




fileType = args.trajectory[-3:]
if fileType == "pdb":
    pdb_trajectory = md.load(args.trajectory, stride=1)
elif fileType == "dcd":
    pdb_trajectory = md.load(args.trajectory, top=input_pdb_filename, stride=1)
    # may use iterload if loading still too slow
else:
    print(f"Unknown fileType {fileType}")
# pdb_trajectory = read_trajectory_pdb_positions(args.trajectory)


# import args.params as params
# default is the params.py under the same folder of trajectory file
if args.params is None:
    location = os.path.dirname(args.trajectory)
    params_location = os.path.join(location, "params.py")
else:
    params_location = args.params
spec = importlib.util.spec_from_file_location("params", params_location)
params = importlib.util.module_from_spec(spec)
spec.loader.exec_module(params)

#Create the merged system
import simtk.openmm
pdb=simtk.openmm.app.PDBFile('clean.pdb')
top=pdb.topology
coord=pdb.positions
forcefield=simtk.openmm.app.ForceField(ffAWSEM.xml,open3SPN2.xml)
s=forcefield.createSystem(top)


dna=open3SPN2.DNA.fromCoarsePDB('clean.pdb', dna_type='B_curved')
with open('protein.seq') as ps:
    protein_sequence_one=ps.readlines()[0]
protein=ffAWSEM.Protein.fromCoarsePDB('clean.pdb',sequence=protein_sequence_one)
dna.periodic=False
protein.periodic=False


#ffAWSEM.copy_parameter_files()
#Clear Forces from the system
keepCMMotionRemover=True
j=0
for i, f in enumerate(s.getForces()):
    if keepCMMotionRemover and i == 0 and f.__class__ == simtk.openmm.CMMotionRemover:
        # print('Kept ', f.__class__)
        j += 1
        continue
    else:
        # print('Removed ', f.__class__)
        s.removeForce(j)
if keepCMMotionRemover == False:
    assert len(s.getForces()) == 0, 'Not all the forces were removed'
else:
    assert len(s.getForces()) <= 1, 'Not all the forces were removed'
forces={}
for i in range(s.getNumForces()):
    force = s.getForce(i)
    force_name="CMMotionRemover"

#Add DNA-protein interaction forces
#AMHgoProteinDNA_flag = False
#US_flag = False

#AMHgoProteinDNA_flag = True
#US_flag = True

for force_name in open3SPN2.protein_dna_forces:
    if force_name in ['AMHgoProteinDNA']:
        force = open3SPN2.protein_dna_forces[force_name](dna, protein, chain_protein='A', chain_DNA='B', k_amhgo_PD=0.8*kilocalorie_per_mole, sigma_sq=0.0225*nanometers**2, aaweight=True)
        s.addForce(force)
        forces.update({force_name: force})
        print(force_name)
    elif force_name in ['String_length_ProteinDNA']:
#        force = open3SPN2.protein_dna_forces[force_name](dna, protein, chain_protein=['A'], chain_DNA=['B'], prot_seg=True, prot_group=[12,15,40,54,56,91,94,120,159,195], DNA_seg=True, DNA_group=list(range(8,13)))
        force = open3SPN2.protein_dna_forces[force_name](dna, protein, chain_protein=['A'], chain_DNA=['B'], prot_seg=True, prot_group=[12,15,40,54,56,91,94,120,159,195])
        s.addForce(force)
        forces.update({force_name: force})
        print(force_name)


# start simulation
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(top, s, integrator, platform)



forceGroupTable = {"Q":1, "Dist":17, "Bond":6, "Angle":7, "Stacking":8, "Dihedral":9, "BasePair":10, "CrossStacking":11, "Exclusion":12, "Electrostatics":13, "ExclProDNA":14, "ElecProDNA":15, "AMHgoProDNA":16, "Con":20, "Rama":21, "Contact":22, "Fragment":23, "Beta":27, "Pap":28, "Helical":29, "DH":30, "Total":list(range(5, 31))}

showEnergy = ["AMHgoProDNA", "Dist"]
# print("Steps", *showEnergy)
if args.output is None:
    outFile = os.path.join(os.path.dirname(args.trajectory), "amh.dat")
else:
    outFile = args.output
with open(outFile, "w") as out:
    line = " ".join(["{0:<8s}".format(i) for i in ["Steps"] + showEnergy])
    print(line)
    out.write(line+"\n")
    # for step, pdb in enumerate(pdb_trajectory):
    #     simulation.context.setPositions(pdb.positions)
    for step in range(len(pdb_trajectory)):
        simulation.context.setPositions(pdb_trajectory.openmm_positions(step))
        e = []
        for term in showEnergy:
            if type(forceGroupTable[term]) == list:
                g = set(forceGroupTable[term])
            elif forceGroupTable[term] == -1:
                g = -1
            else:
                g = {forceGroupTable[term]}
            state = simulation.context.getState(getEnergy=True, groups=g)
            if term == "Q" or term == "Dist":
                termEnergy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
            else:
                termEnergy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
            e.append(termEnergy)
    #     print(*e)
        line = " ".join([f"{step:<8}"] + ["{0:<8.3f}".format(i) for i in e])
        print(line)
        out.write(line+"\n")
    #         print(forceGroupTable[term], state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
