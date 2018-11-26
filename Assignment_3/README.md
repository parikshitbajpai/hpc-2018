# README #

6030G: High Performance Computing

Assignment 2: Langevin Dynamics of Non-Interacting Particles 

Submitted to Prof. Lennaert Van Veen

Submitted by Parikshit Bajpai (Student ID - 100693928)


### Source Codes ###

This repository consists of the following source codes contained in the directory src:

* LangDyn_2D_NonInter.f90 for computing the trajectories and root mean square displacement using a serial code. 

* LangDyn_2D_NonInter_OMP.f90 for computing the trajectories and root mean square displacement using a parallel code implemented using OpenMP directives.

* LangDyn_2D_NonInter_BC_VC.f90 for computing velocity correlation function in addition to the trajectories and root mean square displacement.

### Shell Script ###

* The bash shell script to run all the test cases is named run.sh and can be executed using
```
	  bash ./run.sh
	  ```

### Latex Files ###
* The tex file, bib file and required images are in the folder Latex.
* The report is in the Assignment_2 folder and named Assignment_2_Parikshit.