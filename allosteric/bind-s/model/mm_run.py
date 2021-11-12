#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
from time import sleep
import fileinput
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

from openmmawsem import *
from helperFunctions.myFunctions import *
import ffAWSEM

# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulations")

parser.add_argument("--to", default="./", help="location of movie file")
parser.add_argument("-t", "--thread", type=int, default=-1, help="default is using all that is available")
parser.add_argument("-p", "--platform", type=str, default="OpenCL")
parser.add_argument("--params", type=str, default="params.py")
parser.add_argument("-s", "--steps", type=float, default=1e5, help="step size")
parser.add_argument("-m", "--simulation_mode", type=int, default=0,
                help="default 0: constant temperature,\
                        1: temperature annealing")
args = parser.parse_args()


do = os.system
cd = os.chdir

with open('commandline_args.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


# simulation_platform = "CPU"  # OpenCL, CUDA, CPU, or Reference
# simulation_platform = "OpenCL"
simulation_platform = args.platform
platform = Platform.getPlatformByName(simulation_platform)
if simulation_platform == "CPU":
    if args.thread != -1:
        platform.setPropertyDefaultValue("Threads", str(args.thread))
    print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")




# import args.params as params
spec = importlib.util.spec_from_file_location("params", args.params)
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

protein.k_awsem = 1.0
print("k_awsem:", protein.k_awsem)
dna.k_3spn2 = 1.0
print("k_3spn2:", dna.k_3spn2)

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

##Add 3SPN2 forces
pin_flag = False
#pin_flag = True

for force_name in open3SPN2.forces:
    if force_name in ['PinDNA']:
        if not pin_flag: continue
        force = open3SPN2.forces[force_name](dna, chain_DNA='B', k_pin=400*4.184, pin_file="../../pin_DNA.dat")
        s.addForce(force)
    else:
        force = open3SPN2.forces[force_name](dna)
        if force_name in ['BasePair','CrossStacking']:
            force.addForce(s)
        else:
            s.addForce(force)
    forces.update({force_name: force})
    print(force_name)


#Add AWSEM forces
openAWSEMforces = dict(Connectivity=ffAWSEM.functionTerms.basicTerms.con_term,
                       Chain=ffAWSEM.functionTerms.basicTerms.chain_term,
                       Chi=ffAWSEM.functionTerms.basicTerms.chi_term,
                       Excl=ffAWSEM.functionTerms.basicTerms.excl_term,
                       rama=ffAWSEM.functionTerms.basicTerms.rama_term,
                       rama_pro=ffAWSEM.functionTerms.basicTerms.rama_proline_term,
                       #rama_ss=ffAWSEM.functionTerms.basicTerms.rama_ssweight_term,
                       contact=ffAWSEM.functionTerms.contactTerms.contact_term,
                       beta1 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_1,
                       beta2 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_2,
                       beta3 = ffAWSEM.functionTerms.hydrogenBondTerms.beta_term_3,
                       pap1 = ffAWSEM.functionTerms.hydrogenBondTerms.pap_term_1,
                       pap2 = ffAWSEM.functionTerms.hydrogenBondTerms.pap_term_2,
                       helical = ffAWSEM.functionTerms.hydrogenBondTerms.helical_term,
                       fragment = ffAWSEM.functionTerms.templateTerms.fragment_memory_term,
                       DH = ffAWSEM.functionTerms.debyeHuckelTerms.debye_huckel_term,
                       AMHgo = ffAWSEM.functionTerms.amhgochainTerm.additive_amhgo_chain_term,
                       qBias = ffAWSEM.functionTerms.qBiasTerms.qbias_term,
                      )
protein.setup_virtual_sites(s)
for force_name in openAWSEMforces:
    print(force_name)
    if force_name in ['contact']:
        force = openAWSEMforces[force_name](protein, withExclusion=False, k_contact=params.k_contact, z_dependent=params.z_dependent, inMembrane=params.inMembrane, k_relative_mem=params.k_relative_mem, periodic=params.periodic)
        print(force.getNumExclusions())
        open3SPN2.addNonBondedExclusions(dna,force)
        print(force.getNumExclusions())
    elif force_name in ['fragment']:
        force = openAWSEMforces[force_name](protein, frag_file_list_file="./frags.mem", UseSavedFragTable=True)
    elif force_name in ['AMHgo']:
        force = openAWSEMforces[force_name](protein, "/scratch/xg23/CPEB3/RBD/ct_ZnF/zz.pdb", 200, k_amhgo=4.184*0.3)
    elif force_name in ['DH']:
        force = openAWSEMforces[force_name](protein, coord_ion_file="/scratch/xg23/CPEB3/RBD/ct_ZnF/Zn.index")
    elif force_name in ['qBias']:
        force = openAWSEMforces[force_name](protein, q0=0, reference_pdb_file="/scratch/xg23/CPEB3/RBD/Qbias/freeRRMs.pdb", inter_domain=True, link=100, k_qbias=4.184*10000)
    else:
        force = openAWSEMforces[force_name](protein)
    s.addForce(force)
    forces.update({force_name: force})

#Add DNA-protein interaction forces
AMHgoProteinDNA_flag = False
US_flag = False
length_flag = False
bias_flag = False

AMHgoProteinDNA_flag = True
#US_flag = True
#bias_flag =True

for force_name in open3SPN2.protein_dna_forces:
    if force_name in ['AMHgoProteinDNA']:
        if not AMHgoProteinDNA_flag: continue
        force = open3SPN2.protein_dna_forces[force_name](dna, protein, chain_protein='A', chain_DNA='B', k_amhgo_PD=0.8*kilocalorie_per_mole, sigma_sq=0.0225*nanometers**2, aaweight=True)
    elif force_name in ['StringProteinDNA']:
        if not US_flag: continue
    elif force_name in ['String_length_ProteinDNA']:
        if not length_flag: continue
    elif force_name in ['BiasProteinDNAamhgo']:
        if not bias_flag: continue
    elif force_name in ['ElectrostaticsProteinDNA']:
        force = open3SPN2.protein_dna_forces[force_name](dna, protein, coord_ion_file="/scratch/xg23/CPEB3/RBD/ct_ZnF/Zn.index")
    else:
        force = open3SPN2.protein_dna_forces[force_name](dna, protein)
    s.addForce(force)
    forces.update({force_name: force})
    print(force_name)


# start simulation
print("Simulation Start")
checkpoint_file = "restart"
checkpoint_reporter_frequency = 1000000
reporter_frequency = 5000

# output the native and the structure after minimization
integrator = CustomIntegrator(0.001)
simulation = Simulation(top, s, integrator, platform)
simulation.context.setPositions(coord)  # set the initial positions of the atoms
#simulation.loadState("../cbind.xml")
simulation.reporters.append(PDBReporter(os.path.join(args.to, "native.pdb"), 1))
simulation.step(int(1))
simulation.minimizeEnergy()  # first, minimize the energy to a local minimum to reduce any large forces that might be present
simulation.step(int(1))

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
integrator.setRandomNumberSeed(22028)
simulation = Simulation(top, s, integrator, platform)
simulation.context.setPositions(coord)
#simulation.loadState("../cbind.xml")
# simulation.context.setVelocitiesToTemperature(300*kelvin) # set the initial velocities of the atoms according to the desired starting temperature
simulation.minimizeEnergy() # first, minimize the energy to a local minimum to reduce any large forces that might be present
simulation.reporters.append(StateDataReporter(stdout, reporter_frequency, step=True, potentialEnergy=True, temperature=True)) # output energy and temperature during simulation
simulation.reporters.append(CheckpointReporter(checkpoint_file, checkpoint_reporter_frequency))  # save progress during the simulation
simulation.reporters.append(PDBReporter(os.path.join(args.to, "movie.pdb"), reporter_frequency))  # output PDBs of simulated structures

print("Simulation Starts")
start_time = time.time()

if args.simulation_mode == 0:
    simulation.step(int(args.steps))
elif args.simulation_mode == 1:
    for i in range(100):
        integrator.setTemperature(3*(200-i)*kelvin)
        simulation.step(int(args.steps)/100)


time_taken = time.time() - start_time  # time_taken is in seconds
hours, rest = divmod(time_taken,3600)
minutes, seconds = divmod(rest, 60)
print(f"---{hours} hours {minutes} minutes {seconds} seconds ---")
