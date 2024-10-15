.. dependencies

Dependencies
============

To compile SWIFT, you will need the following libraries:

HDF5
~~~~

Version 1.10.x or higher is required. Input and output files are stored as HDF5
and are compatible with the GADGET-2 specification. A parallel-HDF5 build
and HDF5 >= 1.12.x is recommended when running over MPI.

MPI
~~~
A recent implementation of MPI, such as Open MPI (v3.x or higher), is required,
or any library that implements at least the MPI 3 standard.

Libtool
~~~~~~~
The build system depends on libtool.

FFTW
~~~~
Version 3.3.x or higher is required for periodic gravity.

ParMETIS or METIS
~~~~~~~~~~~~~~~~~
One is required for domain decomposition and load balancing.

GSL
~~~
The GSL 2.x is required for cosmological integration.

In most cases the configuration script will be able to detect the libraries
installed on the system. If that is not the case, the script can be pointed
towards the libraries' location using the following parameters

.. code-block:: bash

  ./configure --with-gsl=<PATH-TO-GSL>

and similar for the other libaries.

Optional Dependencies
=====================

There are also the following *optional* dependencies:

libNUMA
~~~~~~~
libNUMA is used to pin threads.

TCmalloc/Jemalloc/TBBmalloc
~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCmalloc/Jemalloc/TBBmalloc are used for faster memory allocations when available.

Python
~~~~~~
To run the examples, you will need python 3 and some of the standard scientific libraries (numpy, matplotlib).
Some examples make use of the `swiftsimio <https://swiftsimio.readthedocs.io/en/latest/>`_ library,
which is a dedicated and maintained visualisation and analysis library for SWIFT.

GRACKLE
~~~~~~~
GRACKLE cooling is implemented in SWIFT. If you wish to take advantage of it, you will need it installed.

HEALPix C library
~~~~~~~~~~~~~~~~~
This is required for making light cone HEALPix maps. 

CFITSIO
~~~~~~~
This may be required as a dependency of HEALPix.



