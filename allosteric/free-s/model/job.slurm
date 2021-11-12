#!/bin/bash

#SBATCH --job-name=frbd_s
#SBATCH -o trial.out
#SBATCH -e trial.err
#SBATCH --ntasks=1
#SBATCH --partition=commons
#SBATCH --time=14:00:00
#SBATCH --mem=12G
#SBATCH --export=ALL
#SBATCH --gres=gpu:1

python3 mm_run.py free-s -m 0 -s 8e6
