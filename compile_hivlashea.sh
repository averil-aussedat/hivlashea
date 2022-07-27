#!/bin/bash

HIVLASHEA_HOME=.
#Libraries to link for HDF5 : -lsz -lz -lm -lhdf5
#Include path for hdf5.h
#(find your right path by typing "locate hdf5.h" in a terminal, maybe preceded by "updatedb")
INCLUDE_HDF5=/usr/include/hdf5/openmpi
#Library path for libhdf5.a
#(find your right path by typing "locate libhdf5.a" in a terminal, maybe preceded by "updatedb")
LIBRARY_HDF5=/usr/lib/x86_64-linux-gnu/hdf5/openmpi

compile_gcc() {
  export OMPI_CC=gcc
  mpicc -I$HIVLASHEA_HOME/include -I$INCLUDE_HDF5 -L$LIBRARY_HDF5 vp_1d1v_cart_2sp_cemracs2022.c -lsz -lz -lm -lhdf5 -O3 -march=native -std=gnu11 -Wall -Wfloat-conversion -o build/hivlashea.out
}

compile_gcc

