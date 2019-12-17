#!/bin/sh

#SBATCH --job-name=sparta							# job name
#SBATCH -p andrade										# queue
#SBATCH -A jgu-cbdm
#SBATCH --time=01:30:00									# time
#SBATCH --error=./%A.err							# error file name
#SBATCH --output=./%A.out						# output file name
#SBATCH -N 4
#SBATCH -n 256
#SBATCH --mem-per-cpu=1500MB

module load mpi/OpenMPI/3.1.4-GCC-8.3.0

srun --unbuffered -n 256 --mpi=pmi2 python sparta.py ./sc_examples/H3K4me3_hg38_5kb/BC9244609.bed H3K4me3 hg38 5kb RFC-100
