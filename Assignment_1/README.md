# README #

6030G: High Performance Computing

Assignment 1: The matrix-vector product.

Submitted to Prof. Lennaert Van Veen

Submitted by Parikshit Bajpai (Student ID - 100693928)


### Source Codes ###

This repository consists of the following five source codes which are in the directory src:

* matvec.f90 for double loop implementation of matrix vector multiplication. This program runs the double loop for different matrix sizes.

* matvec_fms.f90 for double loop implementation of matrix vector multiplication with matrix size fixed to 16384.

* matvec_OMP.f90 for double loop implementation with parallelisation using OpenMP directives. The matrix size was fixed to 16384.

* matvec_BLAS.f90 for BLAS library implementation of matrix vector multiplication. The matrix size was fixed to 16384.

* matvec_omp_BLAS.f90 for BLAS library implementation with parallelisation using OpenMP directives. The matrix size was fixed to 16384.

* The matrix size can be edited by editing the variable p in the source code. The matrix size would be equal to n = 2^p.


### Shell Script ###

* The bash shell script to run all the test cases is named run.sh and can be executed using
```
	  bash ./run.sh
	  ```
* Note: The shell commands include the paths to the OpenMP directives and  BLAS library and might be different for the machine the code is compiled on. Please edit the path if required.

### Latex Files ###
* The tex file, bib file and required images are in the folder Latex.
* The report is in the Assignment_1 folder and named Assignment_1_Parikshit