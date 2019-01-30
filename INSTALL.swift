
                                Building SWIFT
                                ==============

SWIFT is built from a clean source repository using the commands:

   ./autogen
   ./configure
   make

and from a distribution tarball using:

   ./configure
   make

The compiler choice is GCC by default, but that can be changed using the "CC"
environment variable. This can be just set, or passed on the ./configure
command line, i.e.:

   bash:
      export CC=icc
      ./configure

   [t]csh:
      setenv CC=icc
      ./configure

or:

   ./configure CC=icc

to use an Intel compiler. The main "programs" can be found in the "examples/"
directory. See README for run parameters.

SWIFT has been successfully built and tested with the following compilers:

  - GCC 4.8.x
  - Intel ICC 15.0.x
  - clang 3.4.x

More recent versions and slightly older ones should also be able to
build the software.

It has also been built with Intel and GNU C++ compilers, but that currently
requires the --disable-vec and, for Intel, --disable-compiler-warnings
configure options.

By default an attempt to choose suitable set of optimizing compiler flags
will be made, targeted for the host machine of the build. If this doesn't
work or the binaries will for another architecture then you can stop the
selection of flags using:

   ./configure --disable-optimization

and then supply your own flags using the "CFLAGS" environment variable, as for
CC.

Note that any CFLAGS that you supply will be added to those determined by
configure in all circumstances. To build SWIFT with debugging support you
can use:

    ./configure --enable-debug --disable-optimization

You could also add some additional flags:

    ./configure --enable-debug --disable-optimization CFLAGS="-O2"

for instance. GCC address sanitizer flags can be included using the

    ./configure --enable-sanitizer

option. Note this requires a GCC compiler version of at least 4.8.

By default vectorization is switched on. The highest instruction set
available on the platform will be automatically used. However, not all
implementations of SPH available in the code have vectorized
routines. Vectorization will have to be switched off for these. It can
also be switched off for benchmarking purposes. To do so, you can use:

    ./configure --disable-vec

Please note that to build SWIFT on MacOS, you will need to configure
using

    ./configure --disable-compiler-warnings

due to the incorrect behaviour of the LLVM compiler on this platform
that raises warnings when the pthread flags are passed to the linker.



                                 Dependencies
                                 ============

SWIFT depends on a number of third party libraries that should be available
before you can build it.


 - HDF5:
	A HDF5 library (v. 1.8.x or higher) is required to read and
        write particle data. One of the commands "h5cc" or "h5pcc"
        should be available. If "h5pcc" is located then a parallel
        HDF5 built for the version of MPI located should be
        provided. If the command is not available then it can be
        located using the "--with-hdf5" configure option. The value
        should be the full path to the "h5cc" or "h5pcc" commands.
        SWIFT makes effective use of parallel HDF5 when running on more than
        one node, so this option is highly recommended.

 - MPI:
	To run on more than one node an MPI library that fully
        supports MPI_THREAD_MULTIPLE is required.  Before running configure
        the "mpirun" command should be available in the shell. If your
        command isn't called "mpirun" then define the "MPIRUN"
        environment variable, either in the shell or when running
        configure.

	The MPI compiler can be controlled using the MPICC variable,
	much like the CC one. Use this when your MPI compiler has a
	none-standard name.

 - GSL:
	To use cosmological time integration, a version of the GSL
	must be available.

 - FFTW 3.x:
	To run with periodic gravity forces, a build of the FFTW 3
	library must be available. Note that SWIFT does not make use
	of the parallel capability of FFTW. Calculations are done by
	single MPI nodes independently.

- libtool:
	The build system relies on libtool as well as the other autotools.



                           Optional Dependencies
                           =====================


 - METIS/ParMETIS:
	a build of the METIS or ParMETIS library should be used to
        optimize the load between MPI nodes. This should be found in the
        standard installation directories, or pointed at using the
        "--with-metis" or "--with-parmetis" configuration options.
        In this case the top-level installation directory of the build
        should be given. Note to use METIS or ParMETIS you should supply at
        least "--with-metis". ParMETIS is preferred over METIS when there
        is a choice.

- libNUMA:
	a build of the NUMA library can be used to pin the threads to
        the physical core of the machine SWIFT is running on. This is
        not always necessary as the OS scheduler may do a good job at
        distributing the threads among the different cores on each
        computing node.

        Note that if you have libNUMA outside of the system include
        directories it may fail to compile as the headers do not pass
        the -Wstrict-prototype check of GCC. In that case you will need
        to use --enable-compiler-warnings=yes configure option to stop
        this being an error.

 - tcmalloc / jemalloc / TBBmalloc:
	a build of the tcmalloc library (part of gperftools), jemalloc
	or TBBmalloc can be used be used to obtain faster and more
	scalable allocations than the standard C malloc function part
	of glibc. Using one of these is highly recommended on systems
	with many cores per node. One of the options
	"--with-tcmalloc", "--with-jemalloc" or "--with-tbbmalloc"
	should be passed to the configuration script to use it.

 - gperftools:
	a build of gperftools can be used to obtain good profiling of
        the code. The option "--with-profiler" needs to be passed to
        the configuration script to use it.

 - DOXYGEN:
	the doxygen library is required to create the SWIFT API
        documentation.

 - python:
	Examples and solution script use python and rely on the numpy
	library version 1.8.2 or higher.



                             SWIFT Coding style
                             ==================

The SWIFT source code uses a variation of 'Google' style. The script
'format.sh' in the root directory applies the clang-format-5.0 tool with our
style choices to all the SWIFT C source file. Please apply the formatting
script to the files before submitting a merge request.
