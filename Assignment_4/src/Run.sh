#!/bin/bash
#SBATCH --account=def-veen
#SBATCH --time=0-11:0
#SBATCH --cpus-per-task=32
    
for N_Procs in 1 2 4 8
do
    mpif90 -f90=ifort -heap-arrays -mkl=sequential -o test_1.x main_powermethod.f90
    mpirun -n N_Procs ./test_1.x
done
    



