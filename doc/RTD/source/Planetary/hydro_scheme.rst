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

Note: to access the boundary-improvement method presented in Ruiz-Bonilla+2022,
use the ``planetary_imbalance_RB22`` git branch and compile with
``--with-hydro=planetary-gdf``. However, we instead recommend using the REMIX
SPH scheme, as it has effectively replaced this method.

.. _planetary_remix_hydro:

REMIX SPH
======================

REMIX is an SPH scheme designed to alleviate effects that typically suppress
mixing and instability growth at density discontinuities in SPH simulations
(Sandnes et al. 2025), and also includes the multiple EoS options as the base 
Planetary scheme. For more information on what is included in the REMIX
scheme and how to configure SWIFT to use REMIX, see: :ref:`remix_sph`.
