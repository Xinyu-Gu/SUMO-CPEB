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

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.append(OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

from openmmawsem import *
from helperFunctions.myFunctions import *

# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulations")

parser.add_argument("protein", help="The name of the protein")
parser.add_argument("--name", default="simulation", help="Name of the simulation")
parser.add_argument("--to", default="./", help="location of movie file")
parser.add_argument("-c", "--chain", type=str, default="-1")
parser.add_argument("-t", "--thread", type=int, default=-1, help="default is using all that is available")
parser.add_argument("-p", "--platform", type=str, default="OpenCL")
parser.add_argument("-s", "--steps", type=float, default=1e5, help="step size")
parser.add_argument("-m", "--simulation_mode", type=int, default=0,
                help="default 0: constant temperature,\
                        1: temperature annealing")
parser.add_argument("--params", type=str, default="params.py")
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

proteinName = pdb_id = args.protein
chain=args.chain.upper()
pdb = f"{pdb_id}-openmmawsem.pdb"


if chain == "-1":
    chain = getAllChains(pdb)
    print("Chains to simulate: ", chain)

if args.to != "./":
    # os.system(f"mkdir -p {args.to}")
    os.makedirs(args.to, exist_ok=True)
    os.system(f"cp {args.params} {args.to}/params.py")
    # os.system(f"cp {pdb} {args.to}/{pdb}")
    # pdb = os.path.join(args.to, pdb)

# import args.params as params
spec = importlib.util.spec_from_file_location("params", args.params)
params = importlib.util.module_from_spec(spec)
spec.loader.exec_module(params)


input_pdb_filename = f"{pdb_id}-openmmawsem.pdb"


oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=1.0, chains=chain, seqFromFasta=f"{pdb_id}.fasta", xml_filename=OPENAWSEM_LOCATION+"awsem.xml")  # k_awsem is an overall scaling factor that will affect the relevant temperature scales

# apply forces
forces = [
    con_term(oa),
    chain_term(oa),
    chi_term(oa),
    excl_term(oa, periodic=params.periodic),
    rama_term(oa),
    rama_proline_term(oa),
#    rama_ssweight_term(oa),
    contact_term(oa, withExclusion=False, k_contact=params.k_contact, z_dependent=params.z_dependent, inMembrane=params.inMembrane,
                    k_relative_mem=params.k_relative_mem, periodic=params.periodic),
    beta_term_1(oa),# ,k_beta=4.184*0.5),
    beta_term_2(oa),# ,k_beta=4.184*0.5),
    beta_term_3(oa),# ,k_beta=4.184*0.5),
    pap_term_1(oa),# k_pap=4.184*0.5),
    pap_term_2(oa),# k_pap=4.184*0.5),
    helical_term(oa),# k_helical= 4.184*0.8), 
    fragment_memory_term(oa, frag_file_list_file="./frags.mem", UseSavedFragTable=True),
    debye_huckel_term(oa, coord_ion_file="/scratch/xg23/CPEB3/RBD/ct_ZnF/Zn.index"),
    additive_amhgo_chain_term(oa, "/scratch/xg23/CPEB3/RBD/ct_ZnF/zz.pdb", 200, k_amhgo=4.184*0.3),
    qbias_term(oa, q0=0, reference_pdb_file="/scratch/xg23/CPEB3/RBD/Qbias/freeRRMs.pdb", inter_domain=True, link=100, k_qbias=4.184*10000),
]
oa.addForces(forces)

# start simulation
print("Simulation Start")
checkpoint_file = "restart"
checkpoint_reporter_frequency = 1000000
reporter_frequency = 5000

# output the native and the structure after minimization
integrator = CustomIntegrator(0.001)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
simulation.context.setPositions(oa.pdb.positions)  # set the initial positions of the atoms
#simulation.loadState("relax_restart.xml")
simulation.reporters.append(PDBReporter(os.path.join(args.to, "native.pdb"), 1))
simulation.step(int(1))
simulation.minimizeEnergy()  # first, minimize the energy to a local minimum to reduce any large forces that might be present
simulation.step(int(1))

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
integrator.setRandomNumberSeed(22028)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
simulation.context.setPositions(oa.pdb.positions) # set the initial positions of the atoms
#simulation.loadState("relax_restart.xml")
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
