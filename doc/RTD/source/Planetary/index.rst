.. Planetary Simulations
   Jacob Kegerreis, 14th July 2022

.. _planetary:

Planetary Simulations
=====================

SWIFT is also designed for running planetary simulations
such as of giant impacts, as introduced in
`Kegerreis et al. (2019)  <https://doi.org/10.1093/mnras/stz1606>`_,
and any other non-planetary simulations with more complicated equations
of state and/or multiple materials, etc.

Active development is ongoing of more features and examples,
so please let us know if you are interested in using SWIFT
or have any implementation requests.

You can find an example of creating initial conditions and running an impact
simulation in ``swiftsim/examples/Planetary/`` under ``DemoInitCond/`` and
``DemoImpact/``. Check out their `README.md` files for more info.
The tabulated equation of state files can be downloaded there using
``EoSTables/get_eos_tables.sh``.

To run planetary or similar simulations, SWIFT should be configured to use
the planetary hydrodynamics scheme and equations of state:
``--with-hydro=planetary`` and ``--with-equation-of-state=planetary``.
These allow for multiple materials to be used,
chosen from the several available equations of state.

.. toctree::
   :maxdepth: 1
   :caption: More information:

   Hydro Scheme <hydro_scheme>
   Equations of State <equations_of_state>
   Initial Conditions <initial_conditions>
   Removed Particles <removed_particles>
