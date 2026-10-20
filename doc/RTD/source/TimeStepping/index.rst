.. Details of the time-stepping scheme
   Matthieu Schaller 9th November 2019

.. _time_stepping:

Time integration details
========================

This section of the documentation includes information on the time stepping and
time integration used in SWIFT. It starts with the integer time-line and the
time-bins, then describes how the time-step of a particle is chosen and how the
kick-drift-kick scheme runs as tasks. Then it describes the time-step limiter
and the synchronisation. Finally, the last page explains the content of the
``timesteps.txt`` file.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   integer_time_line
   time_step_criteria
   kick_drift_kick
   timestep_limiter
   timestep_synchronization
   timesteps_file
