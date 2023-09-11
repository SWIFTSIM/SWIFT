.. Running an Example
   Josh Borrow, 5th April 2018
   Mladen Ivkovic, Jan 2023

Running an Example
==================

SWIFT provides a number of examples that you can run in the ``examples/`` directory. 
Many are detailed in their respective ``README`` files, and contain python scripts
(files with the suffix ``.py``) to both generate initial conditions and plot results.
The python scripts usually contain their respective documentation at the top of the 
script file itself.


Sod Shock
~~~~~~~~~

In this example, we will run the 3D SodShock test. You will need to configure and 
compile the code as follows:

.. code-block:: bash
   
   ./configure
   make

Then to run the code, we first download and build the
initial conditions:

.. code-block:: bash

    cd examples/HydroTests/SodShock_3D
    ./getGlass.sh
    python3 makeIC.py
    ../../../swift --hydro --threads=4 sodShock.yml

We can plot the solution with the included python script
as follows: 

.. code-block:: bash

    python3 plotSolution.py 1

The argument ``1`` tells the python plotting script to use the snapshot with number 
1 for the plot.



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
    ../../../swift --cosmology --self-gravity \ 
        --threads=8 small_cosmo_volume_dm.yml


We can plot the solution with the included python script
as follows:

.. code-block:: bash

    python3 plotProjection.py 31

The ``plotProjection.py`` script requires the 
`swiftsimio <https://swiftsimio.readthedocs.io/en/latest/>`_
library.

An example containing both baryonic and dark matter is
``examples/SmallCosmoVolume/SmallCosmoVolume_hydro``. To run with
hydrodynamics, the ``--hydro`` flag needs to be provided as well:

.. code-block:: bash

    ../../../swift --cosmology --self-gravity \ 
        --hydro --threads=8 small_cosmo_volume.yml



