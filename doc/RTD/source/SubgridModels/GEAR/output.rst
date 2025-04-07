.. Sink particles in GEAR model
   Darwin Roduit, 14 July 2024

.. sink_GEAR_model:

Snapshots ouputs
----------------

Here, we provide a summary of the quantities written in the snapshots, in addition to positions, velocities, masses, smoothing lengths and particle IDs.


Sink particles
~~~~~~~~~~~~~~

+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| Name                                  | Description                                 | Units                  | Comments                                          |
+=======================================+=============================================+========================+===================================================+
| ``NumberOfSinkSwallows``              | | Number of merger events with other sinks  | [-]                    |                                                   |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``NumberOfGasSwallows``               | | Number of gas particles accreted          | [-]                    |                                                   |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``TargetMass``                        | | Target mass required to spawn the next    | [U_M]                  | | You can use it to determine if the target mass  |
|                                       | | star particle                             |                        | | is so huge that the sink's mass cannot spawn    |
|                                       |                                             |                        | | such a star. Such rare behaviour may bias the   |
|                                       |                                             |                        | | IMF towards high masses.                        |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``Nstars``                            | | Number of stars created by this sink      | [-]                    |                                                   |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``SwallowedAngularMomentum``          | | Total angular momentum of accreted        | [U_M U_L^2 U_T^{-1}]   |                                                   |
|                                       | | material                                  |                        |                                                   |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``MetalMassFractions``                | | Mass fraction of each tracked metal       | [-]                    | | *Only in GEAR chemistry module.*                |
|                                       | | element                                   |                        | | Array of length ``N`` (number of elements),     |
|                                       |                                             |                        | | set at compile time via                         |
|                                       |                                             |                        | | ``--with-chemistry=GEAR_N``.                    |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``BirthScaleFactors``                 | | Scale factor at the time of sink creation | [-]                    | | Only used in *cosmological* runs.               |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``BirthTimes``                        | | Time when the sink was created            | [U_T]                  | | Only used in *non-cosmological* runs.           |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+



Stars
~~~~~

+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| Name                                  | Description                                 | Units                  | Comments                                          |
+=======================================+=============================================+========================+===================================================+
| ``BirthScaleFactors``                 | | Scale-factors when the stars were born    | [-]                    | | Only used in cosmological runs.                 |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``BirthTimes``                        | | Time when the stars were born             | [U_T]                  | | Only used in non-cosmological runs.             |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``BirthMasses``                       | | Masses of the stars at birth time         | [U_M]                  | | SF and sinks modules                            |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``ProgenitorIDs``                     | | ID of the progenitor sinks or gas         | [-]                    | | SF and sinks modules                            |
|                                       | | particles                                 |                        |                                                   |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``BirthDensities``                    | | Gas density at star formation             | [U_M U_L^{-3}]         | | *Only in SF module*                             |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``BirthTemperatures``                 | | Gas temperature at star formation         | [K]                    | | *Only in SF module*                             |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``Potentials``                        | | Gravitational potential of the star       | [U_L^2 U_T^{-2}]       |                                                   |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``StellarParticleType``               | | Type of the stellar particle:             | [-]                    | | 0: (discrete) single star                       |
|                                       | |                                           |                        | | 1: continuous IMF part star                     |
|                                       | |                                           |                        | | 2: single population star                       |
|                                       | |                                           |                        | | The last type corresponds to legacy IMF stars.  |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``MetalMassFractions``                | | Mass fraction of each metal element       | [-]                    | | *Only in GEAR chemistry module*.                |
|                                       | |                                           |                        | | Array of length ``N`` (number of elements),     |
|                                       | |                                           |                        | | set at compile time by                          |
|                                       | |                                           |                        | | ``--with-chemistry=GEAR_N``.                    |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+

Gas particles
~~~~~~~~~~~~~

Each hydro scheme writes its own set of outputs. Thus, we will only provide the outputs that GEAR writes for gas particles. 
