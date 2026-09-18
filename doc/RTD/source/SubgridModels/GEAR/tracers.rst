.. GEAR sub-grid model
   Darwin Roduit, 10th April 2026

.. _gear_tracers:

Particle tracers
~~~~~~~~~~~~~~~~

Over the course of the simulation, the gas, sink and star particles record some information about their evolution. These are updated for a given particle every time it is active. The GEAR tracers module is located in the directory ``src/tracers/GEAR/``. To enable GEAR tracers, add ``--with-tracers=GEAR`` to your configuration options.

Currently, GEAR tracers are implemented for gas, star, sink and black hole particles.

Gas tracers
-----------

The gas particles record the stellar feedback they receive over their lifetime, separately for supernovae and stellar winds.

+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| Name                                  | Description                                 | Units                       | Comments                                          |
+=======================================+=============================================+=============================+===================================================+
| ``CumulativeMomentumFromSupernovae``  | | Norm of the momentum received from        | [U_M U_L U_T^{-1}]          | | Physical. Scalar sum of the norms, so           |
|                                       | | supernovae, summed over events            |                             | | isotropic kicks do not cancel.                  |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``CumulativeMomentumFromWinds``       | | Norm of the momentum received from        | [U_M U_L U_T^{-1}]          | | Physical. Same convention as for                |
|                                       | | stellar winds, summed over events         |                             | | supernovae.                                     |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``CumulativeEnergyFromSupernovae``    | | Specific internal energy received from    | [U_L^2 U_T^{-2}]            | | Physical.                                       |
|                                       | | supernovae, summed over events            |                             |                                                   |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``CumulativeEnergyFromWinds``         | | Specific internal energy received from    | [U_L^2 U_T^{-2}]            | | Physical. Can be negative if the gas            |
|                                       | | stellar winds, summed over events         |                             | | moved towards the star before the kick.         |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``MaxKickVelocityFromSupernovae``     | | Largest velocity kick received from a     | [U_L U_T^{-1}]              | | Physical.                                       |
|                                       | | single supernova event                    |                             |                                                   |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``MaxKickVelocityFromWinds``          | | Largest velocity kick received from a     | [U_L U_T^{-1}]              | | Physical.                                       |
|                                       | | single stellar winds event                |                             |                                                   |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+

Star tracers
------------

The star particles record their own supernova events, separately for type II and type Ia supernovae. For population particles, an event is a time-step with supernovae, and the number of events can be fractional.

+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| Name                                  | Description                                 | Units                       | Comments                                          |
+=======================================+=============================================+=============================+===================================================+
| ``NumberOfSNIIEvents``                | | Number of SNII events of the star         | [-]                         | | Fractional for population particles.            |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``DensityAtLastSNIIEvent``            | | Gas density around the star at its        | [U_M U_L^{-3}]              | | Physical. 0 if no event yet.                    |
|                                       | | last SNII event                           |                             |                                                   |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``TimeAtLastSNIIEvent``               | | Time of the last SNII event               | [U_T]                       | | Without cosmology. 0 if no event yet.           |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``ScaleFactorAtLastSNIIEvent``        | | Scale factor of the last SNII event       | [-]                         | | With cosmology. 0 if no event yet.              |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``NumberOfSNIaEvents``                | | Number of SNIa events of the star         | [-]                         | | Always 0 for single stars.                      |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``DensityAtLastSNIaEvent``            | | Gas density around the star at its        | [U_M U_L^{-3}]              | | Physical. 0 if no event yet.                    |
|                                       | | last SNIa event                           |                             |                                                   |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``TimeAtLastSNIaEvent``               | | Time of the last SNIa event               | [U_T]                       | | Without cosmology. 0 if no event yet.           |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| ``ScaleFactorAtLastSNIaEvent``        | | Scale factor of the last SNIa event       | [-]                         | | With cosmology. 0 if no event yet.              |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+

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

Black hole tracers
-------------------

The same accretion-rate tracer is available for black hole particles, using ``Snapshots:recording_triggers_bpart`` and the ``num_snapshot_triggers_bpart`` array length (default: 3, matching the sink tracer).

+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
| Name                                  | Description                                 | Units                       | Comments                                          |
+=======================================+=============================================+=============================+===================================================+
| ``AveragedAccretionRates``            | | Accretion rates averaged over the         | [U_M U_T^{-1}]              | | Averaged over the period set by the             |
|                                       | | snapshot trigger intervals                |                             | | first N snapshot triggers.                      |
+---------------------------------------+---------------------------------------------+-----------------------------+---------------------------------------------------+
