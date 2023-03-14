.. Running an Example
   Josh Borrow, 5th April 2018
   Mladen Ivkovic, Jan 2023

Running an Example
==================

Now that you have built the code, you will want to run an example! 
Let's start with a hydrodynamics only example, the 3D Sod Shock.


Sod Shock
~~~~~~~~~

To run the Sod Shock example, you need to follow the following instructions 
(requires ``python3`` with the ``h5py`` and other standard scientific packages, 
as well as ``wget`` for grabbing the glass).


.. code-block:: bash
   
   cd examples/HydroTests/SodShock_3D
   ./getGlass.sh
   python3 makeIC.py
   ../swift --hydro --threads=4 sodShock.yml
   python3 plotSolution.py 1


This will run the 'SodShock' in 3D and produce a nice plot that shows you
how the density has varied. Try running with GIZMO-MFV (this will take
*significantly* longer than with SPH) to see the difference. For that, you
will need to reconfigure with the following options:

.. code-block:: bash
   
   ./configure \
   --with-hydro=gizmo-mfv \
   --with-riemann-solver=hllc

To see the results that you should get, you should check out our developer
wiki at https://gitlab.cosma.dur.ac.uk/swift/swiftsim/wikis/Sod_3D.

If you don't get these results, please contact us on our GitHub page at
https://github.com/SWIFTSIM/swiftsim/issues.




Small Cosmological Volume
~~~~~~~~~~~~~~~~~~~~~~~~~

As a second example, we run a small cosmolgical 
volume containing dark matter only starting at redshift :math:`z = 50`.
Like for the Sod Shock example, it suffices to configure (``./configure``) and 
compile (``make``) the code without any extra flags.

After downloading the initial conditions, we run the code with cosmology and
self-gravity:

.. code-block:: bash

    cd examples/SmallCosmoVolume/SmallCosmoVolume_DM
    ./getIC.sh
    ../../../swift --cosmology --self-gravity --threads=8 small_cosmo_volume_dm.yml


We can plot the solution with the included python script
as follows:

.. code-block:: bash

    python3 plotProjection.py 31


The ``plotProjection.py`` script requires the `swiftsimio <https://swiftsimio.readthedocs.io/en/latest/>`_
library, which is a dedicated and maintained visualisation and analysis
library for SWIFT.


An additional example containing both baryonic and dark matter is
``examples/SmallCosmoVolume/SmallCosmoVolume_hydro``. To run with
hydrodynamics, the ``--hydro`` flag needs to be provided as well:

.. code-block:: bash

    cd examples/SmallCosmoVolume/SmallCosmoVolume_hydro
    ./getIC.sh
    ../../../swift --cosmology --self-gravity --hydro --threads=8 small_cosmo_volume.yml

The solution can again be plotted using the ``plotProjection.py`` script.

