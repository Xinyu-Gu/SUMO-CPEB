#!/bin/bash

#SBATCH --job-name=rrm-f
#SBATCH -o trial.out
#SBATCH -e trial.err
#SBATCH --ntasks=1
#SBATCH --partition=commons
#SBATCH --time=24:00:00
#SBATCH --mem=12G
#SBATCH --export=ALL
#SBATCH --gres=gpu:1


python3 mm_run.py freeRRMs -m 0 -s 6e6
