#!/bin/bash
#SBATCH -A CSC266
#SBATCH -J MiniVATES.jl-bixbyite
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 1:00:00
#SBATCH -p batch-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1

module load PrgEnv-cray-amd julia

export JULIA_DEPOT_PATH=/lustre/polis/csc266/scratch/4pf/julia_depot
cd /lustre/polis/csc266/scratch/4pf/MiniVATES.jl

#export ENABLE_JITPROFILING=1

srun -n $SLURM_NTASKS \
	-c $SLURM_CPUS_PER_TASK \
	--gpus-per-task=$SLURM_GPUS_PER_TASK \
	julia --project test/bixbyite_topaz.jl

	#rocprof --hip-trace --hsa-trace --roctx-trace --sys-trace -o prof.csv \
