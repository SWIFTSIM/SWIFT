#!/bin/bash

#  GNU/GCC toolchain modules with serial HDF5.
#  Source this file.
#
#  Peter W. Draper 20-JAN-2026.

echo "Loading the GNU/GCC toolchain for SWIFT."

module purge
module load gnu_comp/14.1.0
module load openmpi/5.0.3
module load gsl/2.8
module load hdf5/1.14.4
module load parmetis/4.0.3
module load fftw/3.3.10
module load utils
module load llvm/20.1.0
export CLANG_FORMAT_CMD="clang-format"
