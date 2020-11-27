.. Planetary Simulations
   Jacob Kegerreis, 3rd October 2020

.. _planetary:
   
Planetary Simulations
=====================

SWIFT is also designed for running planetary simulations
with a current focus on giant impacts, as introduced in 
`Kegerreis et al. (2019) <https://doi.org/10.1093/mnras/stz1606>`_, MNRAS 487:4.

Active development is ongoing of more features and examples for planetary 
simulations, so please let us know if you are interested in using SWIFT 
or have any implementation requests. 

You can find an example simulation in ``swiftsim/examples/Planetary/``
under ``EarthImpact/``.
The tabulated equation of state files can be downloaded using 
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
