#!/bin/bash
#SBATCH --account=def-veen
#SBATCH --time=0-11:0
#SBATCH --cpus-per-task=32
    

for OMP_NUM_THREADS in 1 2 4 8 16 32
do
    export OMP_NUM_THREADS=$OMP_NUM_THREADS
    ifort -O3 -heap-arrays -qopenmp -o test.x LangDyn_2D_Inter_RBC_DD_OMP.f90
    ./test.x
    cat wtime >> wtime_master_job1
done
    



