#!/bin/bash

#SBATCH --job-name=bs2-frbd
#SBATCH -o trial-re.out
#SBATCH -e trial-re.err
#SBATCH --ntasks=1
#SBATCH --partition=commons
#SBATCH --time=24:00:00
#SBATCH --mem=12G
#SBATCH --export=ALL
#SBATCH --gres=gpu:1

python3 ../relax.py s2-frbd 
