#!/bin/bash

#SBATCH --job-name=5-s2rbd
#SBATCH -o trial.out
#SBATCH -e trial.err
#SBATCH --ntasks=1
#SBATCH --partition=commons
#SBATCH --time=12:00:00
#SBATCH --mem=12G
#SBATCH --export=ALL
#SBATCH --gres=gpu:1

python3 mm_run.py -m 0 -s 8e6
