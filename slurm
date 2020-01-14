#!/bin/sh

#SBATCH --job-name=SIMPA			# job name
#SBATCH -p <PARTITION>				# partition
#SBATCH -A <ACCOUNT>				# the users account
#SBATCH --time=01:30:00				# time
#SBATCH --error=./%A.err			# error file name
#SBATCH --output=./%A.out			# output file name
#SBATCH -N 4					# using 4 nodes
#SBATCH -n 160					# expecting 40 cores per node
#SBATCH --mem-per-cpu=1GB			# memory allocated for job: 160 * 1GB 

# load an installed Open MPI
module load mpi/OpenMPI/3.1.4-GCC-8.3.0

# running the example provided in the repositiory
srun --unbuffered -n 160 --mpi=pmi2 python SIMPA.py --bed ./scExamples/H3K4me3_hg38_5kb/BC8791969.bed --targets H3K4me3

