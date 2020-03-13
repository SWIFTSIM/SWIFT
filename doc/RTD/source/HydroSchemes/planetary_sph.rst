.. Planetary SPH
    Jacob Kegerreis, 13th March 2020

.. _planetary_sph:

Planetary (Density-Energy, Multi-Material) SPH
==============================================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

This scheme is based on the :ref:`minimal` scheme but also allows multiple 
materials, meaning that different SPH particles can be assigned different 
:ref:`equation_of_state` (EoS). 

The Balsara viscosity switch is also used, but can be disabled by 
compiling SWIFT with ``make CFLAGS=-DPLANETARY_SPH_NO_BALSARA``.

To use the planetary scheme and the corresponding planetary EoS, use 

.. code-block:: bash

    ./configure --with-hydro=planetary --with-equation-of-state=planetary

Every SPH particle then requires and carries the additional ``MaterialID`` flag 
from the initial conditions file. This flag indicates the particle's material 
and which EoS it should use. 

See :ref:`planetary` for other related information.
