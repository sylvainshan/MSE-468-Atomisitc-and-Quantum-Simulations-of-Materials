#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks 18
#SBATCH --cpus-per-task 1
#SBATCH --time=00:10:00
#SBATCH --mem=30000

# Bugfix for oneapi issue
# NOTE: This is only needed on helvetios, uncomment if running there.
export FI_PROVIDER=verbs

#These are the libraries necessary for our code to run
# NOTE: This is only needed on helvetios, uncomment them if running there.
module purge
module load intel/2021.6.0  intel-oneapi-mpi/2021.6.0 quantum-espresso/7.2-mpi


PW_LAUCHER="srun" # Use this line for helvetios
# PW_LAUCHER="mpirun" # Use this for the Virtual Machine

# ---------------
# Run
# ---------------

echo "SCF"
$PW_LAUCHER pw.x -input scf.in | tee scf.out

echo "BANDS"
$PW_LAUCHER pw.x -input bands.in | tee bands.out

echo "BANDS.x"
$PW_LAUCHER bands.x -input bandsx.in | tee bandsx.out 

