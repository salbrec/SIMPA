#!/bin/sh

#SBATCH --job-name=SPARTA			# job name
#SBATCH -p <PARTITION>				# partition
#SBATCH -A <ACCOUNT>				# the users account
#SBATCH --time=01:30:00				# time
#SBATCH --error=./%A.err			# error file name
#SBATCH --output=./%A.out			# output file name
#SBATCH -N 4						# using 4 nodes
#SBATCH -n 256						# expecting 64 cores per node
#SBATCH --mem-per-cpu=500MB			# memory allocated for job: 256 * 500MB 

# load an installed Open MPI
module load mpi/OpenMPI/3.1.4-GCC-8.3.0

# running the example provided in the repositiory
srun --unbuffered -n 256 --mpi=pmi2 python sparta.py ./sc_examples/H3K4me3_hg38_5kb/BC9244609.bed H3K4me3 hg38 5kb RFC-100
