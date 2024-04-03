#!/bin/bash
#SBATCH --job-name=term
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:59:00
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --mem=64G


module load julia/1.9.1
#echo 
srun --export=ALL julia -t $NUM_THREADS run.jl $SLURM_ARRAY_TASK_ID "della"
