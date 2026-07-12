.. ShadowSWIFT (Moving mesh hydrodynamics)
   Yolan Uyttenhove September 2023

ShadowSWIFT (moving mesh hydrodynamics)
=======================================

.. warning::
    The moving mesh hydrodynamics solver is currently in the process of being merged into master and will **NOT**
    work on the master branch. To use it, compile the code using the ``moving_mesh`` branch.

This is an implementation of the moving-mesh finite-volume method for hydrodynamics in SWIFT.
To use this scheme, a Riemann solver is also needed. Configure SWIFT as follows:

.. code-block:: bash

    ./configure --with-hydro="shadowswift" --with-riemann-solver="hllc"


Current status
~~~~~~~~~~~~~~

Due to the completely different task structure compared to SPH hydrodynamics, currently only a subset of the features of
SWIFT is supported in this scheme.

-   Hydrodynamics is fully supported in 1D, 2D and 3D and over MPI.

-   Both self-gravity and external potentials are supported.

-   Cosmological time-integration is supported.

-   Cooling and chemistry are supported, with the exception of the ``GEAR_diffusion`` chemistry scheme. Metals are
    properly according to mass fluxes.

-   Choice between periodic, reflective, open, inflow and vacuum boundary conditions (for non-periodic boundary
    conditions, the desired variant must be selected in ``const.h``). Additionally, reflective boundary conditions
    are applied to SWIFT's boundary particles. Configure with ``--with-boundary-particles=<N>`` to use this (e.g. to
    simulate walls).


Caveats
~~~~~~~
These are currently the main limitations of the ShadowSWIFT hydro scheme:

-   Unlike SPH the cells of the moving mesh must form a partition of the entire simulation volume. This means that there
    cannot be empty SWIFT cells and vacuum must be explicitly represented by zero (or negligible) mass particles.
-   Most other subgrid physics, most notably, star formation and stellar feedback are not supported yet.
-   No MHD schemes are supported.
-   No radiative-transfer schemes are supported.
