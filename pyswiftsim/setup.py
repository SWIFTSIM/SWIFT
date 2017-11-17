#!/usr/bin/env python3

__author__ = "loic hausammann"
__copyright__ = "GPLv3"

descr = "Swift is a cosmological hydrodynamical code: SPH With Inter-dependent Fine-grained Tasking"

# import
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy
import os

# read makefile in order to get libraries
makefile = "../src/Makefile"

# save environment
old = {}
def saveEnviron(key, old):
    """ Save environment variable before modifying them
    Key is the environment variable and old is a dictionary
    containing all the modified values.
    """
    old[key] = ""
    if key in os.environ:
        old[key] = os.environ[key]
    return old
        
old = saveEnviron("CC", old)
old = saveEnviron("LDFLAGS", old)


def getValueFromMakefile(f, key):
    """
    Read a file and return the value of a variable.
    f is an open file and value if the requested key.
    The variable should be matching the following patern:
    KEY = VALUE
    """
    test = key + " = "
    ret = None
    for line in f:
        start = line[:len(test)]
        if test == start:
            ret = line[len(test):]

    f.seek(0)
    if ret is None:
        raise Exception("Unable to find key %s" % key)
    else:
        # remove last character "\n"
        return ret[:-1]

with open(makefile, "r") as f:
    hdf5_root = getValueFromMakefile(f, "H5CC")
    cc = getValueFromMakefile(f, "CC")

# python requirement
install_requires = []

# include
include = [
    numpy.get_include(),
    "../src",
    hdf5_root + "/include"
]

# libraries
lib = ["m", "hdf5", "swiftsim"]

# Extension object required by setup
ext = []

# pointer
tmp = Extension("pyswiftsim.pointer",
                ["pyswiftsim/pointer.pyx"],
                include_dirs=include,
                libraries=lib)
tmp = cythonize(tmp)
ext.extend(tmp)

# cooling wrapper
tmp = Extension("pyswiftsim.cooling",
                ["pyswiftsim/cooling.pyx"],
                include_dirs=include,
                libraries=lib)
tmp = cythonize(tmp)
ext.extend(tmp)

# scripts
scripts = []

# data to copy
data = []

# set environment variables, please save the environment variable before modifying them
# path to libswiftsim
ldflags = "LDFLAGS"
if ldflags not in os.environ:
    os.environ[ldflags] = ""
os.environ[ldflags] += " -L" + hdf5_root + "/lib"


# compiler
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


def restoreEnviron(old):
    for key, value in old.items():
        os.environ[key] = value

restoreEnviron(old)
