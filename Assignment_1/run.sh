# Bash shell script to implement the optimization, parallelisation, and standard library codes.
#
# The script contains a nested loop which is split into two parts - a part to compile and generate executable files with the selected optimizations and parallelization, and, a second part to run the executable
#  Step 1: Compile and run the basic double-loop matrix-vector multiplication code without using OpenMP directives. Compiler - GNU (gfortran) & Optimisation - None (-O0).
#          Find the required matrix size and fix the matrix size for further computations.
#
#  Step 2: Compile and run for the different test cases to compare the different schemes.
#     Loop 1: To select the compiler.
#       The cases will be (i) GNU gfortran & (ii) Intel ifort
#
#           Loop 2: To select compiler optimization
#              The cases are (i) No optimization (-O0) & (ii) Aggressive optimization (-O3)
#
#                       Case 1: Compile and run the double-loop code without parallelisation
#                       Case 2: Compile the double-loop code with OpenMP directives
#                       Loop 3: To select the number of threads
#                         The cases are: (i) 1 thread (ii) 2 threads (iii) 4 threads    
#                              Run the double-loop code with OpenMP directives and different number of threads.
#                       End Loop 3
#                       Case 3: Compile and run the BLAS library code without parallelisation
#
#           End Loop 2
#      End Loop 1

echo "----------------------------------------------------------------------------"
echo "----------------------------------------------------------------------------"
echo "6030G HPC                                                 Assignment 1, 2018"
echo "Parikshit Bajpai                                         Student # 100693928"
echo "----------------------------------------------------------------------------"
echo "----------------------------------------------------------------------------"
echo ""
echo "--------------------------Starting Computation------------------------------"
echo ""
echo "--Compiling and running the double-loop code with GNU compiler without optimisation--"
gfortran -O0 -o test.x src/matvec.f90
./test.x
echo ""
echo "Fixed matrix size to 16384 X 16384"
echo ""
echo "---------Compiling and running codes for different approaches to be studied----------"
echo ""
# The following nested loop is for compiling and generating executable files for the different test cases
for comp in gfortran ifort
do
    echo "---------------------------Compiler: $comp -------------------------------- "
    for opti in O0 O3
    do
	echo "--------------------------Optimisation: $opti ------------------------------ "
	echo ""
	echo "Double loop without parallelisation"
	$comp -$opti -o test.x src/matvec_fms.f90
	./test.x
	echo""
	echo "Double loop with parallelisation"
	$comp -$opti -fopenmp -o test.x src/matvec_OMP.f90
	for threads in 1 2 4
	do
	    echo "   Number of threads: $threads"
	    export OMP_NUM_THREADS=$threads
	    ./test.x
	done
	echo ""
	echo "Matrix-vector multiplication using BLAS library"
	$comp -$opti -o test.x src/matvec_BLAS.f90 -L /usr/lib/liblas -lblas # !!! Make sure the BLAS library path is correct 
	./test.x
	echo ""
    done
    echo ""
done
echo ""
echo "Note: The computer was not swapping during computation."
echo ""
echo "-------------------------------------End-------------------------------------"
