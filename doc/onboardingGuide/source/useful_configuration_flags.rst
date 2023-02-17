.. Configuration Options
   Josh Borrow, 5th April 2018

Useful Configuration Flags
==========================

A description of the available options for all flags including the examples
below can be found by using
``./configure  --help``.

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
Several equations of state are made available with this flag.

``--with-cooling=none``
~~~~~~~~~~~~~~~~~~~~~~~
Several cooling implementations (including GRACKLE) are available.

``--with-ext-potential=none``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Many external potentials are available for use with SWIFT. 

.. You can choose
   between them at compile time. Some examples include a central potential, a
   softened central potential, and a sinusoidal potential. You will need to
   configure, for example, the mass in your parameter file at runtime.


