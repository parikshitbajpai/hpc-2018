#! /bin/sh

clear

echo ""
echo "#-----------Langevin Dynamics - Non-Interacting Particles-------------#"
echo ""
echo "Compiler - Intel ifort"
echo ""

echo "Compiling serial code without optimisation"
ifort -O0 -o run.x LangDyn_2D_NonInter_BC_VC.f90 && echo "Compilation OK"
./run.x

