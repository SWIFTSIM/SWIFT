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

Note: the boundary-improvement method presented in Ruiz-Bonilla+2022 can be
accessed on the ``planetary_imbalance_RB22`` git branch.

.. _planetary_remix_hydro:

REMIX SPH
======================

REMIX is an SPH scheme designed to alleviate effects that typically suppress
mixing and instability growth at density discontinuities in SPH simulations
(Sandnes et al. 2025). For more information on what is included in the REMIX
scheme and how to configure SWIFT to use REMIX, see: :ref:`remix_sph`.
