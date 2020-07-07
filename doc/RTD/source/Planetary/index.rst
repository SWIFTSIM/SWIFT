.. Planetary Simulations
   Jacob Kegerreis, 13th March 2020

.. _planetary:
   
Planetary Simulations
=====================

SWIFT is also designed for running planetary simulations
with a current focus on giant impacts, as presented in 
`Kegerreis et al. (2019) <https://doi.org/10.1093/mnras/stz1606>`_, MNRAS 487:4.

More features for planetary simulations are in active development
so please let us know if you are interested in using SWIFT 
or have any implementation requests. For example a new equation of state
or extensions to the tools for creating initial conditions.
    
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
   :maxdepth: 2
   :caption: More information:
   
   Hydro Scheme <hydro_scheme>
   Equations of State <equations_of_state>
   Initial Conditions <initial_conditions>
