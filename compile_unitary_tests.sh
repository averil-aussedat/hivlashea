#!/bin/bash

HIVLASHEA_HOME=.
#Libraries to link for HDF5 : -lsz -lz -lm -lhdf5
#Libraries to link for FFTW : -lfftw3

#Include path for fftw3.h
#(not necessary on local, needed on distant computers)
INCLUDE_FFTW=/
#Library path for -lfftw3
#(not necessary on local, needed on distant computers)
LIBRARY_FFTW=/

#Include path for hdf5.h
#(find your right path by typing "locate hdf5.h" in a terminal, maybe preceded by "updatedb")
INCLUDE_HDF5=/usr/include/hdf5/openmpi
#Library path for libhdf5.a
#(find your right path by typing "locate libhdf5.a" in a terminal, maybe preceded by "updatedb")
LIBRARY_HDF5=/usr/lib/x86_64-linux-gnu/hdf5/openmpi

compile_gcc() {
  export OMPI_CC=gcc
  mpicc -I$HIVLASHEA_HOME/include -I$INCLUDE_HDF5 -L$LIBRARY_HDF5 -I$INCLUDE_FFTW -L$LIBRARY_FFTW unitary_tests.c -lsz -lz -lm -lhdf5 -O3 -march=native -std=gnu11 -Wall -Wfloat-conversion -o build/unitary_tests.out
}

compile_gcc

