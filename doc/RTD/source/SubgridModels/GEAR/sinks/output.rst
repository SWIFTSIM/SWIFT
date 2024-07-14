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
| ``MetalMassFractions``                | | Mass fraction of each metal       | [-]       | | The number of elements ``N`` is determined at   |
|                                       | | element                           |           | | compile time by ``--with-chemistry=GEAR_N``.    |
+---------------------------------------+-------------------------------------+-----------+---------------------------------------------------+


Stars
~~~~~

+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| Name                                  | Description                         | Units     | Comments                    |
+=======================================+=====================================+===========+=============================+
| ``EnergiesReceivedFromJetFeedback``   | | The total energy this particle    | [U_M U_L  |                             |
|                                       | | has received by being kicked as   | ^2 U_t^-2]|                             |
|                                       | | part of jet feedback              |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``LastAGNJetFeedbackScaleFactors``    | | Scale factor when the particle    | [-]       |                             |
|                                       | | was last kicked as part of jet    |           |                             |
|                                       | | feedback, if with cosmology       |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``LastAGNJetFeedbackTimes``           | | Times when the particle was last  | [U_t]     |                             |
|                                       | | last kicked as part of jet        |           |                             |
|                                       | | feedback, if without cosmology    |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``KickedByJetFeedback``               | | How many times this particle has  | [-]       |                             |
|                                       | | been kicked as                    |           |                             |
|                                       | | part of jet feedback              |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
