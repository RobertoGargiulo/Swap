#!/bin/bash
#
gfortran -c -O  exp_sparse.f90 
gfortran -c -O  kyrlov_time_ev.f90
gfortran -c -O  mataid.f90
gfortran -c -O  expokit.f90
gfortran  -O    exp_sparse.o kyrlov_time_ev.o expokit.o mataid.o -o exp_sparse -llapack -lblas ##/lib64/libblas.so.3 
