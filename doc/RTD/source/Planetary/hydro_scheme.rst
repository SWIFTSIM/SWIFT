.. Planetary SPH
    Jacob Kegerreis, 13th March 2020

.. _planetary_hydro:

Planetary Hydro Scheme
======================

This scheme is based on :ref:`minimal` but also allows multiple materials, 
meaning that different SPH particles can be assigned different 
`equations of state <equations_of_state.html>`_ (EoS). 
Every SPH particle then requires and carries the additional ``MaterialID`` flag 
from the initial conditions file. This flag indicates the particle's material 
and which EoS it should use. 

The Balsara viscosity switch is used by default, but can be disabled by 
compiling SWIFT with ``make CFLAGS=-DPLANETARY_SPH_NO_BALSARA``.
