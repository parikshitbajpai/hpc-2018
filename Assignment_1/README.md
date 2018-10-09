# README #

6030G: High Performance Computing
Assignment 1: The matrix-vector product.
Parikshit Bajpai


### Source Codes ###

This repository consists of the following four source codes:

*  matvec.f90 for double loop implementation of matrix vector multiplication. This program runs the double loop for different matrix sizes.

* matvec_omp.f90 for double loop implementation with parallelisation using OpenMP directives. The matrix size was fixed to 8192.

* matvec_BLAS.f90 for BLAS library implementation of matrix vector multiplication. The matrix size was fixed to 8192.

* matvec_omp_BLAS.f90 for BLAS library implementation with parallelisation using OpenMP directives. The matrix size was fixed to 8192.

* The matrix size can be edited by editing the variable p in the source code. The matrix size would be equal to n = 2^p.


### Shell Script ###

* The bash shell script to run all the test cases is named script.sh and can be executed using
```
	  bash ./script.sh
	  ```
* Note: The shell commands include the paths to the OpenMP directives and  BLAS library and might be different for the machine the code is comiled on. Please edit the path if required.

### Latex Files ###
* The tex file, bib file and required images are in the folder Report.
* The report is in the Assignment_1 folder and named Assignment_1_Parikshit