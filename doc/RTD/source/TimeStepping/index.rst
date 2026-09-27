.. Details of the time-stepping scheme
   Matthieu Schaller 9th November 2019

.. _time_stepping:

Time integration details
========================

This section of the documentation includes information on the time
stepping and time integration used in SWIFT. It starts with the integer
time-line and the time-bins, then describes how the time-step of a particle
is chosen and how the kick-drift-kick scheme runs as tasks. The time-step
limiter and the synchronisation are described next. The last page explains the
table of the steps that SWIFT prints, and what the time-step parameters cost.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   integer_time_line
   time_step_criteria
   kick_drift_kick
   timestep_limiter
   timestep_synchronization
   tuning
