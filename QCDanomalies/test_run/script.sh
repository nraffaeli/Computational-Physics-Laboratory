#!/usr/bin/bash
#SBATCH --job-name=ym1_12 # something sensible
#SBATCH --partition=lab4    # this partition is reserved for the lab
#SBATCH --ntasks=16         # the number of MPI processes
#SBATCH --cpus-per-task=1   # this will be the number of OpenMP threads for each tasks / MPI processes
#SBATCH --mem=30G # this is the amount of memory required per node
#SBATCH --hint=nomultithread
#SBATCH --out %x.%j.out     # the standard output of the bash script is redirected to a file
#SBATCH --err %x.%j.err     # and the same for the standard error, see `man sbatch` for more info

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun --cpus-per-task=$SLURM_CPUS_PER_TASK ./ym1_12 -i ym1_12.in -nloc

