.. Configuration Options
   Josh Borrow, 5th April 2018

Configuration Options
=====================

There are many configuration options that SWIFT makes available; a few key
ones are summarised here.

Note that these need to be ran with ``./configure x`` where ``x`` is the
configuration flag.

A description of the available options of the below flags can be found by using
``./configure  --help``.

``--enable-portable-binary``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
By default the compiler flags target the machine ``configure`` runs on, which
can produce a binary that will not run on a different model of node. This
option keeps the tuning but restricts the instruction set, for GCC and clang,
and omits the ``-x`` flag altogether for the Intel compilers. See
:doc:`compiling_code` for when this matters.

``--with-gcc-arch=<arch>``
~~~~~~~~~~~~~~~~~~~~~~~~~~
Name the target architecture explicitly, for example
``--with-gcc-arch=cascadelake``, instead of letting ``configure`` detect the
build machine. Useful when the compute nodes differ from the login node. Read
by the GCC and clang paths only; it has no effect with ``icc`` or ``icx``,
where the equivalent is to pass ``-x...`` in ``CFLAGS``.

``--with-hydro=sphenix``
~~~~~~~~~~~~~~~~~~~~~~~~
There are several hydrodynamical schemes available in SWIFT. You can choose
between them at compile-time with this option.

``--with-riemann-solver=none``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Some hydrodynamical schemes, for example GIZMO, require a Riemann solver.

``--with-kernel=cubic-spline``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Several kernels are made available for use with the hydrodynamical schemes.
Choose between them with this compile-time flag.

``--with-hydro-dimension=3``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Run problems in 1, 2, and 3 (default) dimensions.

``--with-equation-of-state=ideal-gas``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Several equations of state are made available with this flag. Also consider
``--with-adiabatic-index``.

``--with-cooling=none``
~~~~~~~~~~~~~~~~~~~~~~~
Several cooling implementations (including GRACKLE) are available.

``--with-ext-potential=none``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Many external potentials are available for use with SWIFT. You can choose
between them at compile time. Some examples include a central potential, a
softened central potential, and a sinusoidal potential. You will need to
configure, for example, the mass in your parameter file at runtime.


