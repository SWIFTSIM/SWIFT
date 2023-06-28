.. Compiling the Code
   Josh Borrow, 5th April 2018


Compiling SWIFT
===============

Compilers
---------

SWIFT is a C99 code, and as such requires a C compiler that is able
to work with code built for that standard.

SWIFT has been tested with the Intel, GCC, LLVM (clang) C compilers.

We suggest:

+ Intel >= 2018
+ GCC >= 8.2.0
+ LLVM >= 7.0.0

We have specific issues with the following compilers:

+ GCC 7.3.0 with the ``-mskylake-avx512`` flag.

Dependencies
------------

To compile SWIFT, you will need the following libraries:

HDF5
~~~~

Version 1.8.x or higher is required. Input and output files are stored as HDF5
and are compatible with the existing GADGET-2 specification. Please consider
using a build of parallel-HDF5, as SWIFT can leverage this when writing and
reading snapshots. We recommend using HDF5 > 1.10.x as this is *vastly superior*
in parallel.

HDF5 is widely available through system package managers.

MPI
~~~
A recent implementation of MPI, such as Open MPI (v2.x or higher), is required,
or any library that implements at least the MPI 3 standard.
MPI implementations are widely available through system package managers.

Running SWIFT on OmniPath atchitechtures with Open MPI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When running SWIFT on an OmniPath system we suggest that Open MPI v3.1.3 or higher
is used. A bug in the ``psm2`` library causes communications to be lost. It is
possible to run SWIFT with older versions (tested with v2.1.x) of Open MPI so
long as ``psm`` is used instead of ``psm2``, i.e. that you invoke ``mpirun``
with ``--mca btl vader,self -mca mtl psm``.

Libtool
~~~~~~~
The build system depends on libtool. Libtool is widely available through system 
package managers.

FFTW
~~~~
Version 3.3.x or higher is required for periodic gravity. FFTW  is widely available
through system package managers or on http://fftw.org/.

ParMETIS or METIS
~~~~~~~~~~~~~~~~~
One of these libraries is required for domain decomposition and load balancing. 
Source codes for them libraries are available 
`here for METIS <https://github.com/KarypisLab/METIS>`_ and 
`here for ParMETIS <https://github.com/KarypisLab/ParMETIS>`_ .

GSL
~~~
The GSL is required for cosmological integration. GSL is widely available through
system package managers.


Optional Dependencies
---------------------

There are also the following *optional* dependencies.

libNUMA
~~~~~~~
libNUMA is used to pin threads (but see INSTALL.swift).

TCmalloc/Jemalloc
~~~~~~~~~~~~~~~~~
TCmalloc/Jemalloc are used for faster memory allocations when available.

DOXYGEN
~~~~~~~
You can build documentation for SWIFT with DOXYGEN.

Python
~~~~~~
To run the examples, you will need python 3 and some of the standard scientific 
libraries (numpy, matplotlib). Some examples make use of the 
`swiftsimio <https://swiftsimio.readthedocs.io/en/latest/>`_ library.

GRACKLE
~~~~~~~
GRACKLE cooling is implemented in SWIFT. If you wish to take advantage of it, you 
will need it installed. It can be found `here <https://github.com/grackle-project/grackle>`_.


HEALPix C library
~~~~~~~~~~~~~~~~~~~

This is required for making light cone HEALPix maps. Note that by default HEALPix 
builds a static library which cannot be used to build the SWIFT shared library. 
Either HEALPix must be built as a shared library or -fPIC must be added to the C 
compiler flags when HEALPix is being configured.

CFITSIO
~~~~~~~

This may be required as a dependency of HEALPix.


Initial Setup
-------------

We use autotools for setup. To get a basic running version of the code use:

.. code-block:: bash

  ./autogen.sh
  ./configure
  make

the executable binaries are found in the top directory.

MacOS Specific Oddities
~~~~~~~~~~~~~~~~~~~~~~~

To build on MacOS you will need to disable compiler warnings due to an
incomplete implementation of pthread barriers. DOXYGEN also has some issues on
MacOS, so it is best to leave it out. To configure:

.. code-block:: bash

  ./configure --disable-compiler-warnings --disable-doxygen-doc

When using the clang compiler, the hand-written vectorized routines
have to be disabled. This is done at configuration time by adding
the flag ``--disable-hand-vec``.

Trouble Finding Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~

If the configure script is having trouble finding your libraries for you, it
may be that they are in nonstandard locations. You can link the specific
library locations by using ``--with-<LIBRARY>=<PATH>``. For example for the
HDF5 library,

.. code-block:: bash
   
   ./configure --with-hdf5=/path/to/hdf5_root

More information about what needs to be provided to these flags is given in
``./configure --help``.
