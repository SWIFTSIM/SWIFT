.. Planetary SPH
    Jacob Kegerreis, 3rd February 2019

.. _planetary_sph:

Planetary (Density-Energy, Multi-Material) SPH
==============================================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

This scheme is the same as the Minimal SPH scheme but also allows multiple 
materials, meaning that different SPH particles can be assigned different 
:ref:`equation_of_state` (EoS).

To use the planetary scheme and the corresponding planetary EoS, use 

.. code-block:: bash

    ./configure --with-hydro=planetary --with-equation-of-state=planetary

Every SPH particle then requires and carries the additional ``MaterialID`` flag 
from the initial conditions file. This flag indicates the particle's material 
and which EoS it should use. 