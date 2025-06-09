Initial Setup
=============

We use autotools for setup. To get a basic running version of the code (the executable binaries are found in the top directory), use:

.. code-block:: bash

  ./autogen.sh
  ./configure
  make


MacOS Specific Oddities
~~~~~~~~~~~~~~~~~~~~~~~

To build on MacOS you will need to enable compiler warnings due to an
incomplete implementation of pthread barriers. DOXYGEN also has some issues on
MacOS, so it is best to leave it out. To configure:

.. code-block:: bash

  ./configure --enable-compiler-warnings \
      --disable-doxygen-doc

When using the ``clang`` compiler, the hand-written vectorized routines
have to be disabled. This is done at configuration time by adding
the flag ``--disable-hand-vec``.


