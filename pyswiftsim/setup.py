#!/usr/bin/env python3

__author__ = "loic hausammann"
__copyright__ = "GPLv3"

descr = "Swift is a cosmological hydrodynamical code: SPH With Inter-dependent Fine-grained Tasking"

# import
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy
import os

makefile = "../src/Makefile"

def getValueFromMakefile(f, value):
    test = value + " = "
    ret = None
    for line in f:
        start = line[:len(test)]
        if test == start:
            ret = line[len(test):]

    f.seek(0)
    if ret is None:
        raise Exception("Unable to find value %s" % value)
    else:
        # remove last character "\n"
        return ret[:-1]

# python requirement
install_requires = []

# include
with open(makefile, "r") as f:
    # need to remove \n
    hdf5_include = getValueFromMakefile(f, "H5CC") + "/include"
    cc = getValueFromMakefile(f, "CC")
    
include = [
    numpy.get_include(),
    "../src",
    hdf5_include
]

# libraries
lib = ["m", "hdf5", "-L../src/.libs/", "swiftsim"]

# Extension object required by setup
ext = []


# cooling wrapper
tmp = Extension("cooling",
                ["src/cooling.pyx"],
                include_dirs=include,
                libraries=lib)
tmp = cythonize(tmp)
ext.extend(tmp)

# scripts
scripts = []

# data to copy
data = ["../src/.libs/libswiftsim.so"] 
old_cc = ""
if "CC" in os.environ:
    old_cc = os.environ['CC']

os.environ["CC"] = cc

setup(
    name = "pyswiftsim",
    version      = "0.6.0",
    author       = "Hausammann Loic",
    author_email = "loic.hausammann@epfl.ch",
    description  = descr,
    license      = "GPLv3",
    keywords     = "swift hpc cosmology",
    url          = "https://gitlab.cosma.dur.ac.uk/swift/swiftsim",

    packages         = find_packages(),

    data_files       = data,

    scripts          = scripts,

    install_requires = install_requires,

    ext_modules      = ext,

)

os.environ["CC"] = old_cc
