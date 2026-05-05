.. GEAR sub-grid model
   Darwin Roduit, 10th April 2026

.. _gear_tracers:

Particle tracers
~~~~~~~~~~~~~~~~

Over the course of the simulation, the gas, sink and star particles record some information about their evolution. These are updated for a given particle every time it is active. The GEAR tracers module is located in the directory ``src/tracers/GEAR/``. To enable GEAR tracers, add ``--with-tracers=GEAR`` to your configuration options.

Currently, GEAR tracers are only implemented for sink particles.

Sink tracers
------------

The tracers track the long-term evolution of the sinks' accretion rates and SFR between different time slices. These tracers are used to capture the high-frequency variability of the the accretion rate and SFR that might be lost between snapshots.

The time slices are specified by ``Snapshots:recording_triggers_sink``. By default, the number of recording time slices is 3. You can change this value in ``src/tracers_triggers.h`` and then recompile the code. These outputs are arrays of length ``num_snapshot_triggers_sink``.

+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| Name                                  | Description                                 | Units                       | Comments                                          |
+=======================================+=============================================+=============================+===================================================+
| ``AveragedAccretionRates``            | | Accretion rates averaged over the         | [U_M U_T^{-1}]              | | Averaged over the period set by the             |
|                                       | | snapshot trigger intervals                |                             | | first N snapshot triggers.                      |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``AveragedStarFormationRates``        | | Star formation rates averaged over the    | [U_M U_T^{-1}]              | | Averaged over the period set by the             |
|                                       | | snapshot trigger intervals                |                             | | first N snapshot triggers.                      |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
