#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64

source /public3/soft/modules/module.sh
module load mpi/oneAPI/2022.1
module load mpi/intel/2022.1

srun main
