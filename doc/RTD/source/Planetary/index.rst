.. Planetary Simulations
   Jacob Kegerreis, 14th July 2022

.. _planetary:

Planetary Simulations
=====================

SWIFT is also designed for running planetary simulations
such as of giant impacts, as introduced in Kegerreis+2019,
and any other types of simulations with more complicated equations of state
and/or multiple materials, etc.

Active development is ongoing of more features and examples,
so please let us know if you are interested in using SWIFT
or have any implementation requests.

You can find an example simulation in ``swiftsim/examples/Planetary/``
under ``EarthImpact/``, as well as some hydro tests for comparison with other
schemes. The tabulated equation of state files can be downloaded using
``EoSTables/get_eos_tables.sh``.

Planetary simulations are currently intended to be run with
SWIFT configured to use the planetary hydrodynamics scheme and equations of state:
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
