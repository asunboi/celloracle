#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=72000
#SBATCH --time=40:00:00

cd $SLURM_SUBMIT_DIR

python songpeng.py
