#!/bin/bash

#  Intel toolchain modules with serial HDF5 and metis.
#  Source this file.

module purge
module load intel_comp/2024.2.0 compiler-rt tbb compiler mpi
module load hdf5/1.14.4
module load fftw/3.3.10
module load gsl/2.8
module load metis/4.0.3
module load python/3.12.4
module load llvm/20.1.0
export CLANG_FORMAT_CMD="clang-format"
