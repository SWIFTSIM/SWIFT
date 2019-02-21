.. Compiling the Code
   Josh Borrow, 5th April 2018


Compiling SWIFT
===============

Dependencies
------------

To compile SWIFT, you will need the following libraries:

HDF5
~~~~

Version 1.8.x or higher is required. Input and output files are stored as HDF5
and are compatible with the existing GADGET-2 specification. Please consider
using a build of parallel-HDF5, as SWIFT can leverage this when writing and
reading snapshots. We recommend using HDF5 > 1.10.x as this is `vastly superior`
in parallel.

MPI
~~~
A recent implementation of MPI, such as Open MPI (v2.x or higher), is required,
or any library that implements at least the MPI 3 standard.

Running SWIFT on OmniPath atchitechtures with Open MPI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When running SWIFT on an OmniPath system we suggest that Open MPI v3.1.3 or higher
is used. A bug in the ``psm2`` library causes communications to be lost. It is
possible to run SWIFT with older versions (tested with v2.1.x) of Open MPI so
long as ``psm`` is used instead of ``psm2``, i.e. that you invoke ``mpirun``
with ``--mca btl vader,self -mca mtl psm``.

Libtool
~~~~~~~
The build system depends on libtool.

FFTW
~~~~
Version 3.3.x or higher is required for periodic gravity.

ParMETIS or METIS
~~~~~~~~~~~~~~~~~
One is required for domain decomposition and load balancing.

libNUMA
~~~~~~~
libNUMA is used to pin threads (but see INSTALL.swift).

GSL
~~~
The GSL is required for cosmological integration.


Optional Dependencies
---------------------

There are also the following _optional_ dependencies.

TCmalloc/Jemalloc
~~~~~~~~~~~~~~~~~
TCmalloc/Jemalloc are used for faster memory allocations when available.

DOXYGEN
~~~~~~~
You can build documentation for SWIFT with DOXYGEN.

Python
~~~~~~
To run the examples, you will need python and some of the standard scientific libraries (numpy, matplotlib). Some examples use Python 2 scripts, but the more recent ones use Python 3 (this is specified in individual READMEs).

GRACKLE
~~~~~~~
GRACKLE cooling is implemented in SWIFT. If you wish to take advantage of it, you will need it installed.


Initial Setup
-------------

We use autotools for setup. To get a basic running version of the code
(the binary is created in swiftsim/examples) on most platforms, run

.. code-block:: bash

  ./autogen.sh
  ./configure
  make


MacOS Specific Oddities
~~~~~~~~~~~~~~~~~~~~~~~

To build on MacOS you will need to disable compiler warnings due to an
incomplete implementation of pthread barriers. DOXYGEN also has some issues on
MacOS, so it is best to leave it out. To configure:

.. code-block:: bash

  ./configure --disable-compiler-warnings --disable-doxygen-doc


Trouble Finding Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~

If the configure script is having trouble finding your libraries for you, it
may be that they are in nonstandard locations. You can link the specific
library locations by using ``--with-<LIBRARY>=<PATH>``. For example for the
HDF5 library,

.. code-block:: bash
   
   ./configure --with-hdf5=/path/to/h5cc

More information about what needs to be provided to these flags is given in
``./configure --help``.
