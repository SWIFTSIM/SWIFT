.. Sink particles in GEAR model
   Darwin Roduit, 14 July 2024

.. sink_GEAR_model:

.. _gear_output_fields:

Snapshots ouputs
----------------

Here, we provide a summary of the quantities written in the snapshots, in addition to positions, velocities, masses, smoothing lengths and particle IDs.

Most tracer outputs are summarised on the :ref:`gear_tracers` page. The gas and star fields the radiation model registers are grouped by channel below instead: radiation pressure's own diagnostics are tracers-module fields; HII photoionization's and the interstellar radiation field's are not.

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
| ``BirthDensities``                    | | Physical densities at the time of birth   | [U_M U_L^{-3}]         | | Stored at birth time/redshift.                  |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``BirthTemperatures``                 | | Temperatures at the time of birth         | [U_K]                  | | Stored at birth time/redshift.                  |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``AccretionRates``                    | | Physical instantaneous accretion rates    | [U_M U_T^{-1}]         | | At the current step                             |
+---------------------------------------+---------------------------------------------+------------------------+---------------------------------------------------+
| ``StarFormationRates``                | | Star formation rates of the particles     | [U_M U_T^{-1}]         | | If negative, stores last time/scale-factor      |
|                                       |                                             |                        | | at which gas was star-forming.                  |
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

Since hydro scheme writes its own set of outputs, we only provide the outputs that ``GEAR`` writes for gas particles. 

+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| Name                                     | Description                                                 | Units                  | Comments                                           |
+==========================================+=============================================================+========================+====================================================+
| ``SmoothedMetalMassFractions``           | | Mass fraction of each metal element                       | [-]                    | | *Only in GEAR chemistry module.*                 |
|                                          | | smoothed over the SPH kernel                              |                        | | Array of length ``N``, set at compile time by    |
|                                          |                                                             |                        | | ``--with-chemistry=GEAR_N``                      |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``MetalMassFractions``                   | | Raw (non-smoothed) mass fraction of                       | [-]                    | | *Only in GEAR chemistry module.*                 |
|                                          | | each metal element                                        |                        | | Same layout as above.                            |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``HI``                                   | | Mass fraction of neutral H (:math:`\mathrm{H}`)           | [-]                    | | *Only if* ``GRACKLE_1 to 3``                     |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``HII``                                  | | Mass fraction of ionized H (:math:`\mathrm{H}^+`)         | [-]                    | | *Only if* ``GRACKLE_1 to 3``                     |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``HeI``                                  | | Mass fraction of neutral He (:math:`\mathrm{He}`)         | [-]                    | | *Only if* ``GRACKLE_1 to 3``                     |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``HeII``                                 | | Mass fraction of singly ionized He (:math:`\mathrm{He}^+`)| [-]                    | | *Only if* ``GRACKLE_1 to 3``                     |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``HeIII``                                | | Mass fraction of doubly ionized He                        | [-]                    | | *Only if* ``GRACKLE_1 to 3``                     |
|                                          | | (:math:`\mathrm{He}^{++}`)                                |                        |                                                    |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``e``                                    | | Free electron mass fraction (:math:`\mathrm{e}^-`)        | [-]                    | | *Only if* ``GRACKLE_1 to 3``                     |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``HM``                                   | | Mass fraction of :math:`\mathrm{H}^-`                     | [-]                    | | *Only if* ``GRACKLE_2 or 3``                     |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``H2I``                                  | | Mass fraction of neutral :math:`\mathrm{H}_2`             | [-]                    | | *Only if* ``GRACKLE_2 or 3``                     |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``H2II``                                 | | Mass fraction of ionized :math:`\mathrm{H}_2^+`           | [-]                    | | *Only if* ``GRACKLE_2 or 3``                     |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``DI``                                   | | Mass fraction of neutral D (:math:`\mathrm{D}`)           | [-]                    | | *Only if* ``GRACKLE_3``                          |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``DII``                                  | | Mass fraction of ionized D (:math:`\mathrm{D}^+`)         | [-]                    | | *Only if* ``GRACKLE_3``                          |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+
| ``HDI``                                  | | Mass fraction of :math:`\mathrm{HD}`                      | [-]                    | | *Only if* ``GRACKLE_3``                          |
+------------------------------------------+-------------------------------------------------------------+------------------------+----------------------------------------------------+

Black holes
~~~~~~~~~~~

See :ref:`gear_black_holes` for the physics behind these fields. The averaged accretion-rate tracer is documented separately, on the :ref:`gear_tracers` page.

.. list-table::
   :header-rows: 1
   :widths: 30 45 15 30

   * - Name
     - Description
     - Units
     - Comments
   * - ``DynamicalMasses``
     - Dynamical (gravitational, resolved) mass of the black hole particle
     - [U_M]
     - Generic ``Masses`` is not written for black holes; use this field.
   * - ``SubgridMasses``
     - Subgrid (unresolved, physical) mass of the black hole
     - [U_M]
     - Starts at the dynamical mass at formation and diverges from it as the BH accretes; see :ref:`gear_black_holes`.
   * - ``FormationScaleFactors`` / ``FormationTimes``
     - Scale-factor / time at which the BH was formed
     - [-] / [U_T]
     - Scale-factor in cosmological runs, time otherwise.
   * - ``GasDensities``
     - Co-moving density of the gas around the black hole
     - [U_M U_L^{-3}]
     -
   * - ``GasSoundSpeeds``
     - Co-moving sound-speed of the gas around the black hole
     - [U_L U_T^{-1}]
     -
   * - ``EnergyReservoirs``
     - Energy currently stored in the AGN feedback reservoir
     - [U_M U_L^2 U_T^{-2}]
     -
   * - ``AccretionRates``
     - Physical instantaneous accretion rate
     - [U_M U_T^{-1}]
     -
   * - ``TotalAccretedMasses``
     - Total mass accreted since birth (main progenitor only)
     - [U_M]
     - Excludes mass accreted by any merged-in black holes.
   * - ``CumulativeNumberOfSeeds``
     - Number of BH seeds merged into this black hole (including itself)
     - [-]
     -
   * - ``NumberOfMergers``
     - Number of BH-BH mergers experienced
     - [-]
     - Excludes mergers accumulated by any merged-in black holes.
   * - ``LastHighEddingtonFractionScaleFactors`` / ``LastHighEddingtonFractionTimes``
     - Scale-factor / time the BH last exceeded ``eddington_fraction_for_recording``
     - [-] / [U_T]
     - -1 if never reached.
   * - ``LastMinorMergerScaleFactors`` / ``LastMinorMergerTimes``
     - Scale-factor / time of the last minor merger
     - [-] / [U_T]
     -
   * - ``LastMajorMergerScaleFactors`` / ``LastMajorMergerTimes``
     - Scale-factor / time of the last major merger
     - [-] / [U_T]
     -
   * - ``SwallowedAngularMomenta``
     - Angular momentum accumulated by swallowing gas particles
     - [U_M U_L^2 U_T^{-1}]
     -
   * - ``AccretedAngularMomenta``
     - Angular momentum accumulated through subgrid accretion
     - [U_M U_L^2 U_T^{-1}]
     -
   * - ``GasRelativeVelocities``
     - Peculiar velocity of the gas relative to the black hole
     - [U_L U_T^{-1}]
     -
   * - ``GasCircularVelocities``
     - Circular velocity of the gas at the smoothing radius
     - [U_L U_T^{-1}]
     -
   * - ``NumberOfSwallows`` / ``NumberOfDirectSwallows``
     - Number of gas particles swallowed, including / excluding those swallowed by merged-in black holes
     - [-]
     - ``NumberOfSwallows`` includes swallows by merged-in black holes.
   * - ``NumberOfRepositions`` / ``NumberOfRepositionAttempts``
     - Number of repositioning events / of steps with an eligible reposition target
     - [-]
     - See :ref:`gear_black_holes` for how to disable repositioning entirely.
   * - ``NumberOfTimeSteps``
     - Number of time steps the black hole was active
     - [-]
     -
   * - ``ViscosityFactors``
     - Suppression factor from the Rosas-Guevara et al. (2015) model
     - [-]
     - Only meaningful if ``with_angmom_limiter`` is enabled.
   * - ``BirthGasDensities``
     - Physical gas density at the time of birth
     - [U_M U_L^{-3}]
     -
   * - ``NumberOfGasNeighbours``
     - Number of gas neighbours within the black hole's kernel
     - [-]
     -
   * - ``FeedbackDeltaT``
     - Temperature increase applied in the most recent AGN feedback event
     - [U_K]
     -
   * - ``LastRepositionVelocities``
     - Speed of the most recent repositioning jump
     - [U_L U_T^{-1}]
     - 0 if never repositioned, or run without a prescribed repositioning speed.
   * - ``NumberOfHeatingEvents``
     - Number of (thermal) energy injections
     - [-]
     - A single AGN event can heat several particles, incrementing this by more than 1.
   * - ``NumberOfAGNEvents``
     - Number of time steps in which the BH did AGN feedback
     - [-]
     -
   * - ``LastAGNFeedbackScaleFactors`` / ``LastAGNFeedbackTimes``
     - Scale-factor / time of the last AGN feedback event
     - [-] / [U_T]
     -
   * - ``AccretionLimitedTimeSteps``
     - Accretion-limited time-step
     - [U_T]
     - The particle's actual time-step may be shorter due to the minimum allowed value.
   * - ``AGNTotalInjectedEnergies``
     - Cumulative energy injected into gas by AGN feedback
     - [U_M U_L^2 U_T^{-2}]
     -
   * - ``AccretionBoostFactors``
     - Booth & Schaye (2009) accretion boost factor
     - [-]
     - Only meaningful if ``with_boost_factor`` is enabled.
   * - ``GasTemperatures``
     - Temperature of the gas surrounding the black hole
     - [U_K]
     -
   * - ``EnergyReservoirThresholds``
     - Minimum reservoir energy required to do feedback
     - [-]
     - In units of the (constant) target heating temperature increase.
   * - ``EddingtonFractions``
     - Accretion rate in units of the Eddington rate
     - [-]
     - Based on the unlimited accretion rate, so can exceed ``max_eddington_fraction``.
   * - ``Potentials``
     - Gravitational potential of the black hole
     - [U_L^2 U_T^{-2}]
     -

.. _gear_output_hii:

HII photoionization
~~~~~~~~~~~~~~~~~~~~

The fields of this and the next two sections exist only in a build that compiles their channel (``./configure --with-subgrid-radiation``, see :ref:`gear_radiation`); ``EffectiveTemperatures`` exists whenever any channel is compiled.

See :ref:`gear_radiation_hii` for the physics. None of the six fields below are tracers-module fields, so none need ``--with-tracers=GEAR``:

.. list-table::
   :header-rows: 1
   :widths: 25 45 15 30

   * - Name
     - Description
     - Units
     - Comments
   * - ``HIIRegionRadii``
     - Comoving radius the star's HII region reached at its last budget rebuild
     - [U_L]
     - Star. The search algorithm's own bookkeeping, not a direct measurement of the gas's thermal state.
   * - ``HIIRegionMasses``
     - Gas mass the star currently holds ionized
     - [U_M]
     - Star. Same caveat as ``HIIRegionRadii``.
   * - ``IsIonizedFlags``
     - Is this gas particle currently flagged as ionized?
     - [-]
     - Gas.
   * - ``HIIStarIDs``
     - ID of the star that ionized this particle
     - [-]
     - Gas.
   * - ``FinalHIIRegionRadii``
     - ``HIIRegionRadii`` retired at the star's death or ineligibility
     - [U_L]
     - Star.
   * - ``FinalHIIRegionMasses``
     - ``HIIRegionMasses`` retired the same way
     - [U_M]
     - Star.

.. _gear_output_radiation_pressure:

Radiation pressure
~~~~~~~~~~~~~~~~~~~~

See :ref:`gear_radiation_pressure` for the physics. These gas fields need ``--with-tracers=GEAR``:

.. list-table::
   :header-rows: 1
   :widths: 25 45 15 30

   * - Name
     - Description
     - Units
     - Comments
   * - ``CumulativeMomentumFromRadiationPressure``
     - Norm of the momentum received from radiation pressure, summed over events (scalar sum of the norms, so isotropic kicks do not cancel)
     - [U_M U_L U_T^{-1}]
     - Gas.
   * - ``MaxKickVelocityFromRadiationPressure``
     - Largest single-event velocity kick this particle received from radiation pressure
     - [U_L U_T^{-1}]
     - Gas.

.. _gear_output_isrf:

Interstellar radiation field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See :ref:`gear_isrf` for the physics. In a build that compiles the ISRF, ``PELuminosities`` and ``LWLuminosities`` are written for every GEAR run with feedback, whether or not the ISRF module or ``--with-tracers=GEAR`` are used; both read 0 unless ``with_interstellar_radiation_field`` is on.

.. list-table::
   :header-rows: 1
   :widths: 25 45 15 30

   * - Name
     - Description
     - Units
     - Comments
   * - ``PELuminosities``
     - Star's current non-ionizing PE-band (6-11.2 eV) luminosity, physical
     - [U_M U_L^2 U_T^{-3}]
     - Star.
   * - ``LWLuminosities``
     - Star's current Lyman-Werner-band (11.2-13.6 eV) luminosity, physical
     - [U_M U_L^2 U_T^{-3}]
     - Star.
   * - ``EffectiveTemperatures``
     - Photospheric effective temperature of the star, from the radiation table. For a population particle, the value at the upper mass bound of the stars still alive (the hottest surviving star), not an IMF average
     - [U_K]
     - Star. Diagnostic only: no feedback channel uses it. 0 with no radiation table, or a table with no ``Teff`` dataset.

The remaining fields are gas fields registered by the GEAR feedback module's own I/O code, not the tracers module, so they do not need ``--with-tracers=GEAR``. All are physical quantities with no scale-factor exponent of their own.

.. list-table::
   :header-rows: 1
   :widths: 25 45 15 30

   * - Name
     - Description
     - Units
     - Comments
   * - ``PESpecificEnergies``
     - Local specific PE-band field
     - [U_L^2 U_T^{-2}]
     -
   * - ``LWSpecificEnergies``
     - Local specific Lyman-Werner-band field
     - [U_L^2 U_T^{-2}]
     -
   * - ``LWPhotonSpecificEnergies``
     - Local Lyman-Werner-band photon-number moment, energy-equivalent at a fixed reference photon energy, not a photon count
     - [U_L^2 U_T^{-2}]
     - Diagnostic only: no heating or dissociation rate reads it.

.. _gear_output_isrf_packet:

With ``ISRF_propagation: 1`` and a build configured with ``--enable-isrf-debugging-checks`` (see :ref:`gear_radiation`) the transport writes more fields, all zero otherwise. A build without that option writes none of the ``Packet*`` fields. Every two-component field holds the PE band, then the Lyman-Werner band. "Ledger" fields close the energy budget of a run: they are the only way to check energy conservation from a snapshot. "Diagnostic" fields are for debugging the transport and are not needed for science. Switch fields off with the output selection (see :ref:`Output_selection_label`); the names carry the particle type as a suffix, for example ``PacketSlotDirections_Gas``. With 12 slots, the diagnostic and duplicate gas fields take 464 bytes per gas particle and the ledger gas fields 144 bytes, against 24 bytes for the three fields above (all double precision). The slot fields (432 bytes) scale with the number of slots.

.. list-table::
   :header-rows: 1
   :widths: 30 40 15 15

   * - Name
     - Description
     - Units
     - Kind
   * - ``PacketKernelMeanSpecificEnergies``, ``PacketKernelMeanSpecificPhotonNumbers``
     - Kernel mean of the own field of the packet. Same values as ``PESpecificEnergies`` and ``LWSpecificEnergies``
     - [U_L^2 U_T^{-2}]
     - Diagnostic (duplicate). Gas.
   * - ``PacketSpecificEnergies``, ``PacketSpecificPhotonNumbers``
     - Own field of the particle before the kernel mean. Two components
     - [U_L^2 U_T^{-2}]
     - Diagnostic. Gas.
   * - ``PacketSlotDirections``
     - Unit direction of the packet of each reception slot (three numbers per slot)
     - none
     - Diagnostic. Gas.
   * - ``PacketSlotOriginDistances``
     - Comoving distance from the origin of the packet to the particle, per slot
     - [U_L]
     - Diagnostic. Gas.
   * - ``PacketSlotInboxEnergies``, ``PacketSlotInboxDirections``
     - Energy and direction in the inbox of each slot
     - [U_M U_L^2 U_T^{-2}], none
     - Diagnostic. Gas.
   * - ``PacketInboxEnergies``, ``PacketInboxPhotonNumbers``
     - Energy received since the last step of the particle
     - [U_M U_L^2 U_T^{-2}]
     - Ledger. Gas.
   * - ``PacketReservoirEnergies``, ``PacketReservoirPhotonNumbers``
     - Energy that waits beyond the arrival bins of the inbox
     - [U_M U_L^2 U_T^{-2}]
     - Ledger. Gas.
   * - ``PacketHeldEnergies``, ``PacketHeldPhotonNumbers``
     - Energy held for sending this step; zero between steps
     - [U_M U_L^2 U_T^{-2}]
     - Ledger. Gas.
   * - ``PacketAbsorbedEnergies``, ``PacketEscapedEnergies``, ``PacketRedshiftedEnergies`` and the ``...PhotonNumbers`` twins
     - Cumulative energy absorbed by dust, left with no downstream neighbour, and removed by the cosmological expansion
     - [U_M U_L^2 U_T^{-2}]
     - Ledger. Gas.
   * - ``PacketEmittedEnergies``, ``PacketEscapedEnergies``, ``PacketRemovedEnergies``, ``PacketInheritedEnergies`` and the ``...PhotonNumbers`` twins
     - Cumulative energy emitted, put in empty emission cones, removed with a swallowed gas particle and inherited, per star
     - [U_M U_L^2 U_T^{-2}]
     - Ledger. Star.
   * - ``PacketRemovedEnergies``, ``PacketInheritedEnergies`` and the ``...PhotonNumbers`` twins
     - As for stars
     - [U_M U_L^2 U_T^{-2}]
     - Ledger. Sink.

The photon-number fields are in energy units: the photon number times the lower edge energy of the band (6.0 eV for PE, 11.2 eV for Lyman-Werner).

Gas that leaves a box with open faces is deleted together with the light it holds, so no particle field keeps its ledger. The run total of that light is written instead to the header of the snapshot, in the group ``SubgridScheme``, by a build with ``--enable-isrf-debugging-checks``. Each of the four attributes has two components (PE band, then Lyman-Werner band), is summed over the MPI ranks and is in snapshot units. ``PacketLeftBoxEscapedEnergies`` and ``PacketLeftBoxEscapedPhotonNumbers`` hold the light in flight in the gas when it left, which the energy balance counts as escaped. ``PacketLeftBoxInheritedEnergies`` and ``PacketLeftBoxInheritedPhotonNumbers`` hold the absorbed, escaped and redshifted light that the departed particles had already booked. Add them to the sums of the particle ledger fields to close the budget of a run with open faces.

Three optional gas input fields let an initial-conditions file seed the local field directly; see :ref:`gear_isrf` for their use.
