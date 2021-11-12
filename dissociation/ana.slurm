#!/bin/bash

#SBATCH --job-name=rrm-b
#SBATCH -o trial-ana.out
#SBATCH -e trial.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=scavenge
#SBATCH --time=00:20:00
#SBATCH --mem=12G
#SBATCH --export=ALL


python3 /scratch/xg23/CPEB3/RBD/ct_ZnF/0.8amhgo/mm_amh.py -p CPU --thread 1
