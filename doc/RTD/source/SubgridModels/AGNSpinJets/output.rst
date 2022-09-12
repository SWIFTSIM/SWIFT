.. AGN spin and jet model
   Filip Husko, 1 April 2022

.. AGN_spin_jet:

Black holes
-----------


+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| Name                                  | Description                         | Units     | Comments                    |
+=======================================+=====================================+===========+=============================+
| ``AngularMomentumDirections``         | | The direction of the angular      | [-]       | | Array of length           |
|                                       | | momentum (spin) vector of the BH  |           | | 3 for each particle       |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``AccretionDiscAspectRatios``         | | Aspect ratio, H/R, of the subgrid | [-]       |                             |
|                                       | | accretion disk surrounding each   |           |                             |
|                                       | | black hole                        |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``AccretionModes``                    | | Type of subgrid accretion disk    | [-]       |                             |
|                                       | | surrounding the black holes       |           |                             |
|                                       | | 0 - Thick disk, 1 - Thin disk,    |           |                             |
|                                       | | 2 - Slim disk                     |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``CosAccretionDiskAngle``             | | Cosine of the angle between the   | [-]       |                             |
|                                       | | spin vector and the gas angular   |           |                             |
|                                       | | momentum vector around the BH     |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``InjectedJetEnergies``               | | Total jet energy injected into    | [U_M U_L  |                             |
|                                       | | surroundings of this BH           | ^2 U_t^-2]|                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``JetEfficiencies``                   | | The efficiency of jet launching,  | [-]       |                             |
|                                       | | i.e. the jet power divided by the |           |                             |
|                                       | | accretion rate times c*c          |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``JetReservoirs``                     | | The remaining jet energy left to  | [U_M U_L  |                             |
|                                       | | be launched from the BH           | ^2 U_t^-2]|                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``JetTimeSteps``                      | | Jet-limited time steps of the BHs | [U_t]     |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``LastAGNJetScaleFactors``            | | Last AGN jet scale factors when   | [-]       |                             |
|                                       | | the BH did jet feedback, if       |           |                             |
|                                       | | cosmology is turned on            |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``LastAGNJetTimes``                   | | Last AGN jet times when the BH    | [U_t]     |                             |
|                                       | | did jet feedback, if cosmology is |           |                             |
|                                       | | turned off                        |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``NumberOfAGNJetEvents``              | | Total number of times this BH did | [-]       |                             |
|                                       | | jet feedback                      |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``NumberOfJetParticlesLaunched``      | | Total number of times this BH     | [-]       |                             |
|                                       | | launched particles as part of jet |           |                             |
|                                       | | feedback                          |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``Spins``                             | | Dimensionless spin parameters of  | [-]       |                             |
|                                       | | the black holes. Negative values  |           |                             |
|                                       | | indicate retrograde accretion     |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+
| ``RadiativeEfficiencies``             | | The efficiency of radiative       | [-]       |                             |
|                                       | | feedback, i.e. the radiative      |           |                             |
|                                       | | power divided by the accretion    |           |                             |
|                                       | | rate times c*c                    |           |                             |
+---------------------------------------+-------------------------------------+-----------+-----------------------------+


Tracers (gas and stars)
-----------------------

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
