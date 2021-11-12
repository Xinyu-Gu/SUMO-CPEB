#!/bin/bash

#SBATCH --job-name=bs2-frbd
#SBATCH -o trial.out
#SBATCH -e trial.err
#SBATCH --ntasks=1
#SBATCH --partition=commons
#SBATCH --time=12:00:00
#SBATCH --mem=12G
#SBATCH --export=ALL
#SBATCH --gres=gpu:1

python3 mm_run.py s2-frbd 
