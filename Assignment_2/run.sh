#-------------------Angevin Dynamics of Non-Interacting Particles-----------------------#
#                                                                                        #
# This script runs the source codes to reproduce the results presented in the report.    #
# The script compiles the serial and parallel codes using Intel's ifort compiler and      #
# finally runs the modified script for computing the velocity correlation.               #
#                                                                                        #
#----------------------------------------------------------------------------------------#


#! /bin/sh

clear

echo ""
echo "#--------------Langevin Dynamics - Non-Interacting Particles----------------#"
echo ""
echo "Compiler - Intel ifort"
echo ""

# Compiling and running the series and parallel codes for comparing the wall times and calculating speed-up and efficiency. 
cd bin
echo "Compiling serial code without optimisation"
ifort -O0 -o run_s.x ../src/LangDyn_2D_NonInter.f90 && echo "Compilation OK"
./run_s.x
echo ""
echo "Compiling parallel code without optimisation"
ifort -O0 -fopenmp -o run_p.x ../src/LangDyn_2D_NonInter_OMP.f90 && echo "Compilation OK"
./run_p.x
echo ""
echo "Compiling serial code with aggressive optimisation"
ifort -O3 -o run_s.x ../src/LangDyn_2D_NonInter.f90 && echo "Compilation OK"
./run_s.x
echo ""
echo "Compiling parallel code with aggressive optimisation"
ifort -O3 -fopenmp -o run_p.x ../src/LangDyn_2D_NonInter_OMP.f90 && echo "Compilation OK"
./run_p.x

# Cleaning up the bin folder
rm *.*

# Compiling and running the modified code for finding the velocity correlation function. 
echo ""
echo "Velocity Correlation - Compiling serial code with aggressive optimisation"
ifort -O3 -o run.x ../src/LangDyn_2D_NonInter_BC_VC.f90 && echo "Compilation OK"
./run_p.x

echo ""
gnuplot Plot

echo ""
echo "----------------------------------------------------------------------------"
echo ""

