.. Sink particles in GEAR model
   Darwin Roduit, 14 July 2024

.. sink_GEAR_model:

Snapshots ouputs
----------------

Here, we provide a summary of the quantities written in the snapshots, in addition to positions, velocities, masses and particle IDs.

Sink particles
~~~~~~~~~~~~~~

+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| Name                                  | Description                         | Units     | Comments                                          |
+=======================================+=====================================+===========+===================================================+
| ``NumberOfSinkSwallows``              | | Number of sink merger events      | [-]       |                                                   |
|                                       | |                                   |           |                                                   |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| ``NumberOfGasSwallows``               | | Number of gas swallowed           | [-]       |                                                   |
|                                       | |                                   |           |                                                   |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| ``TargetMass``                        | | Sink target mass to spawn the     | [U_M]     | | You can use it to determine if the target mass  |
|                                       | | next star particle                |           | | is so huge that the sink's mass cannot spawn    |
|                                       | |                                   |           | | such a star. Such rare behaviour may bias the   |
|                                       | |                                   |           | | IMF towards high masses.                        |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| ``Nstars``                            | | Number of stars created by the    | [-]       |                                                   |
|                                       | | the sink particles                |           |                                                   |
|                                       | |                                   |           |                                                   |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| ``SwallowedAngularMomentum``          | | Total angular momentum swallowed  | [U_M U_L  |                                                   |
|                                       | | by the sink particles             |  ^2 U_T   |                                                   |
|                                       | |                                   |  ^-1]     |                                                   |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| ``MetalMassFractions``                | | Mass fraction of each metal       | [-]       | | Array of length ``N`` for each particles. The   |
|                                       | | element                           |           | | number of elements ``N`` is determined at       |
|                                       | |                                   |           | | compile time by ``--with-chemistry=GEAR_N``.    |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+


Stars
~~~~~

+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| Name                                  | Description                         | Units     | Comments                                          |
+=======================================+=====================================+===========+===================================================+
| ``BirthScaleFactors``                 | | Scale-factors when the stars were | [-]       | Only used in cosmological runs.                   |
|                                       | | born                              |           |                                                   |
|                                       | |                                   |           |                                                   |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| ``BirthTimes``                        | | Time when the stars were          | [U_T]     | Only used in non-cosmological runs.               |
|                                       | | born                              |           |                                                   |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| ``BirthMasses``                       | | Masses of the stars at brith time | [U_M]     |                                                   |
|                                       | |                                   |           |                                                   |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| ``ProgenitorIDs``                     | | ID of the progenitor sinks        | [-]       |                                                   |
|                                       | |                                   |           |                                                   |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
| ``StellarParticleType``               | | Type of the stellar particle:     | [-]       | | The last type correspond to star particles in   |
|                                       | | 0: (discrete) single star         |           | | the previous model, i.e. representing the full  |
|                                       | | 1: continuous IMF part star       |           | | IMF                                             |
|                                       | | 2: single population star         |           | |                                                 |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+
