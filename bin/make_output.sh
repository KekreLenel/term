#!/bin/bash
#SBATCH --job-name=term_output
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:59:00
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --mem=64G


module load julia/1.8.2

srun --export=ALL julia run.jl "output"  "della"
