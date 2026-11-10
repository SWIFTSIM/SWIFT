.. Running an Example
   Josh Borrow, 5th April 2018
   Mladen Ivkovic, Jan 2023

Running an Example
==================

Now that you have built the code, you will want to run an example! 
Let's start with a hydrodynamics only example, the 3D Sod Shock.


Python requirements
~~~~~~~~~~~~~~~~~~~

The examples use python scripts to make their initial conditions and to plot
their results. These need ``python3`` (version 3.11 or later) and the packages
listed in ``examples/requirements.txt``, which can be installed with:

.. code-block:: bash

   python3 -m pip install -r examples/requirements.txt

The main one is `swiftsimio <https://swiftsimio.readthedocs.io/en/latest/>`_
(version 12.4.0 or later), which is used to write the initial conditions and
to read the snapshots. The download scripts also need ``wget``.

A few examples need additional packages:

==================  ========================================================
Package             Examples
==================  ========================================================
``woma`` (< 1.5)    ``Planetary/DemoImpact``, ``Planetary/DemoImpactInitCond``,
                    ``Planetary/JupiterLikePlanet``,
                    ``Planetary/RayleighTaylor_EarthLike_3D``,
                    ``Planetary/RayleighTaylor_JupiterLike_3D``
``p_tqdm``          the ``makeMovie.py`` scripts of ``HydroTests/BlobTest_3D``,
                    ``HydroTests/Rayleigh-Taylor_2D`` and
                    ``HydroTests/SquareTest_2D``
``pyyaml``          ``RadiativeTransferTests/CosmoCoolingTest``,
                    ``RadiativeTransferTests/CosmoHeatingTest``,
                    ``RadiativeTransferTests/HeatingTest``
                    (``generate_output_times.py``)
``yt``              ``HydroTests/KeplerianRing`` (optional plots)
``galpy``           ``GravityTests/NFW_Halo``
``agama``           ``GravityTests/ExternalAgama``
``pNbody``          ``GravityTests/ExternalAgama``,
                    ``IsolatedGalaxy/IsolatedGalaxy_multi_component``
``hmf``             ``SmallCosmoVolume/SmallCosmoVolume_VELOCIraptor``
``healpy``          ``SmallCosmoVolume/SmallCosmoVolume_lightcone``
``velociraptor``    ``SmallCosmoVolume/SmallCosmoVolume_DM`` (CSDS analysis,
                    which also needs SWIFT configured with
                    ``--enable-csds --enable-python``)
==================  ========================================================

Versions 1.5 and later of ``woma`` fail to compute the equations of state used
by these examples, so use an earlier version
(``python3 -m pip install "woma<1.5"``).


Sod Shock
~~~~~~~~~

To run the Sod Shock example, you need to follow the following instructions
(requires the python packages listed above, as well as ``wget`` for grabbing
the glass).


.. code-block:: bash
   
   cd examples/HydroTests/SodShock_3D
   ./getGlass.sh
   python3 makeIC.py
   ../../../swift --hydro --threads=4 sodShock.yml
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

