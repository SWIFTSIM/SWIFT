
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
directory.

By default an attempt to choose suitable set of optimizing compiler flags
will be made, targetted for the host machine of the build. If this doesn't
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

Dependencies: needs to be filled in...


    



