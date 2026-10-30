.. GEAR sub-grid model interstellar radiation field

.. _gear_isrf:

Interstellar radiation field
============================

Young stars produce a local, non-ionizing ultraviolet radiation field that heats the gas through the photoelectric effect on dust grains and dissociates molecular hydrogen. The GEAR interstellar radiation field (ISRF) module follows this field in two bands:

- **PE**, 6 to 11.2 eV, which drives photoelectric heating on dust,
- **LW** (Lyman-Werner), 11.2 to 13.6 eV, which drives :math:`\mathrm{H}_2` photodissociation.

Each star carries a band luminosity interpolated from the radiation datasets of its yields table (``GEARFeedback:yields_table``), as a function of its mass, metallicity and age. The luminosity is deposited on the gas particles inside the star's SPH kernel. Every receiving gas particle then attenuates the field it sees by a dust extinction factor :math:`\exp(-\kappa_\mathrm{eff} \Sigma)`, with the column :math:`\Sigma = \rho \, \ell` built from the particle's own density and a path :math:`\ell` set by the receiver-side extinction mechanism (see below; by default, the star-to-particle separation). The dust opacity scales with the particle's metallicity, so a metal-free gas particle is not shielded.

Injection alone illuminates only the stars' immediate neighbourhoods. Optionally (``GEARFeedback:ISRF_propagation: 1``), the light of each star is instead handed from particle to particle by a packet-hop transport: the light travels as packets of energy from one gas particle to the next, reaches gas beyond the source kernels at a finite speed and is shadowed by the gas in between. This is the only transport scheme of the module. The terms of the transport are explained in :ref:`gear_isrf_glossary`.

The resulting per-particle field is handed to Grackle every step. The Lyman-Werner band alone sets the :math:`\mathrm{H}_2` photodissociation rate (Grackle mode 2 or 3 only; see :ref:`gear_radiation`'s Grackle cooling mode paragraph). Photoelectric heating uses a different field, the Habing field :math:`G_0`. The Habing band spans 6 to 13.6 eV, so :math:`G_0` sums the PE and Lyman-Werner band energies. Lyman-Werner energy therefore still feeds photoelectric heating at every Grackle mode, including below mode 2.

Working configurations are shipped in ``examples/SubgridTests/StellarFeedback/ISRF/``; each example directory carries its own README with its configure line, run command and check script.

Configuring and compiling
-------------------------

The module is part of the GEAR feedback model and is coupled to Grackle, so a run needs at least:

.. code:: bash

   ./configure --with-feedback=GEAR --with-chemistry=GEAR_10 --with-stars=GEAR \
               --with-cooling=grackle_2 --with-grackle=$GRACKLE_ROOT

For the packet-hop transport (see below), add the HEALPix library (``--with-chealpix=<path>``, which can need CFITSIO; see :doc:`the compilation instructions </GettingStarted/compiling_code>`) and, optionally, the number of reception slots:

.. code:: bash

   ./configure --with-feedback=GEAR --with-chemistry=GEAR_10 --with-stars=GEAR \
               --with-cooling=grackle_2 --with-grackle=$GRACKLE_ROOT \
               --with-chealpix=$CHEALPIX_ROOT --with-isrf-packet-slots=12

To keep the energy ledgers of the transport and write them to the snapshot, add ``--enable-isrf-debugging-checks`` (off by default; it makes the particles larger and is independent of ``--enable-debugging-checks``; a restart file can only be resumed by a build with the same setting). Without it the snapshot has no ``Packet*`` field.

See :ref:`gear_radiation` for the per-band Grackle mode gating shared by every radiation channel (and why ``--with-subgrid=GEAR`` is the wrong shortcut here), and :ref:`gear_grackle_cooling` for the cooling module itself.

**The yields table must carry the ISRF datasets.** This module additionally needs the ``L_PE``, ``L_LW``, ``Integrated_L_PE`` and ``Integrated_L_LW`` datasets in the table's ``Data/Radiation`` group, on top of the datasets every radiation channel needs. See :ref:`gear_radiation_tables` for the full requirement and how to get and check a table.

Switching the module on
-----------------------

A single parameter enables the module:

.. code:: YAML

   GEARFeedback:
     with_interstellar_radiation_field: 1   # Master switch of the ISRF module

This switch turns on both channels and forces the Grackle flags they need (the ISRF field, the dust chemistry, the photoelectric heating and, from ``grackle_2`` upwards, the radiative-transfer rate channel) on internally, so you do not set those yourself.

With ``with_interstellar_radiation_field: 1`` and everything else left at its default, stars illuminate their own kernels, the receiving gas is shielded by its own dust column, and no transport takes place. This is enough for the photoelectric heating channel in a well-resolved interstellar medium, and it is the cheapest configuration.

Propagation
-----------

Set ``GEARFeedback:ISRF_propagation: 1`` to transport the light of the stars with the packet hop. The default, ``0``, adds the dose of a star to the gas inside its kernel at once.

A star splits its band luminosity into emission cones and hands each cone to the gas particles it reaches. Each gas particle gathers the light that reaches it in reception slots, which are cells of equal solid angle fixed in a rotated frame. Once per own time step the particle absorbs a part of the light in its dust and hands the rest to its downstream neighbours with a transmission cone :math:`W \max(0, \mu)^p`. The field that Grackle reads is the kernel mean of the light that the particle holds. The transport rides the time step of the gas: there is no separate radiation clock.

The transport also runs in cosmological boxes, see :ref:`gear_isrf_cosmology`.

The transport needs a yields table with the datasets ``MeanPhotonEnergyPE`` and ``Integrated_MeanPhotonEnergyPE`` in the ``Data/Radiation`` group (see :ref:`gear_radiation_tables`). SWIFT stops at start-up if one is missing. Without dust the light of gas at zero metallicity is never absorbed, so check the initial metallicity of the gas.

Light held by a gas particle with no downstream neighbour escapes and leaves the transport. A star counts the light that its emission cones send to no gas as escaped too. Light held by a gas particle that becomes a star, or that a sink swallows, leaves with it and is booked to that star or sink as removed. ``./configure --enable-isrf-debugging-checks`` keeps these energy ledgers and writes them to the snapshot (see :ref:`gear_radiation`).

``GEARFeedback:ISRF_packet_light_speed`` (default ``cap``)
  Speed limit of the light. ``cap``: a share arrives no earlier than the separation divided by the limit speed, in the frame of the gas. ``fixed``: as ``cap``, and the time step of each gas particle is bounded by :math:`f_{\Delta t}\,a\,h / v`, so that the light moves at the limit speed. ``unlimited``: a share arrives at the next step of its receiver, and the light can move faster than :math:`c` (debugging only, SWIFT warns at start-up).

``GEARFeedback:ISRF_packet_light_speed_fraction_of_c`` (default ``0.99``)
  Limit speed :math:`v` as a fraction of :math:`c`, in :math:`(0, 1]`. Read by ``cap`` and ``fixed``. Raise it only to approach :math:`c`; lower it to slow the light on purpose.

``GEARFeedback:ISRF_packet_light_speed_time_step_margin`` (default ``0.5``)
  Margin :math:`f_{\Delta t}` of the time step bound of ``fixed``, in :math:`(0, 1.69]`. Above about :math:`0.8` the bulk of the light lags behind the limit speed while its front keeps it.

``GEARFeedback:ISRF_packet_cone_power`` (default ``4``)
  Power :math:`p \geq 1` of the transmission cone. A larger value gives a narrower beam and less lateral spread.

``GEARFeedback:ISRF_packet_read_mean_time_Myr`` (default ``0``)
  Time constant of a running mean of the field that a particle reads. ``0`` switches it off. Use a few typical time steps of the gas, converted to Myr, to smooth step-to-step noise of the field.

``GEARFeedback:ISRF_packet_octant_frame_angles_deg`` (default ``[24, 36, 66]``) and ``ISRF_packet_octant_frame_seed_offset`` (default ``0``)
  Rotation of the frame of the reception slots, in degrees about x, y and z, and the seed of an extra random rotation (``0``: none). The defaults avoid alignment of the slot boundaries with the axes of a box. Change them only to test that a result does not depend on the frame.

Test beams
~~~~~~~~~~

Two optional parameters inject a pulse of light without a star. They exist to validate the transport and are not meant for a production run.

``GEARFeedback:ISRF_packet_test_beam_1`` and ``GEARFeedback:ISRF_packet_test_beam_2`` (default: off)
  Each is a list of 11 numbers in internal units: the lower corner of a box (3), the upper corner of the box (3), the direction of the beam (3), and the energy given to each gas particle in the PE band and in the Lyman-Werner band (2). Every gas particle inside the box receives that energy once, at its first step, as light that travels along the direction (it is normalised). The pulse then propagates like any other light. A beam with zero energy in both bands is off. The energy is an energy, not a specific energy, in the internal energy unit of the run (internal mass unit times the square of the internal velocity unit). For example, with the units of the example parameter file (:math:`10^{10}\,\mathrm{M}_\odot`, 1 kpc, 1 km/s), one internal energy unit is :math:`1.988\times10^{53}` erg, so a pulse of :math:`10^{45}` erg in the PE band is written as :math:`10^{45} / 1.988\times10^{53} = 5.03\times10^{-9}`.

The injected energy is in no ledger, so a test adds it up from the parameters. SWIFT prints one message for each active beam, stops if an energy is negative, the direction is zero or the upper corner is not above the lower one, and stops if a beam is set with ``ISRF_propagation: 0``. In a cosmological run the corners are comoving coordinates.

.. _gear_isrf_cosmology:

Cosmological runs
~~~~~~~~~~~~~~~~~

The transport needs no extra parameter in a cosmological run: it reads the scale factor :math:`a` and the Hubble rate :math:`H` from the cosmology. A static box has :math:`H = 0`, so none of what follows applies to it.

Units.
  Positions, separations and smoothing lengths of the transport are comoving. The dust absorption, the time steps, the speed of light and the field that Grackle reads are physical. The optical depth of a hop is the dust absorption rate times the physical path, :math:`a` times the comoving one. The field handed to Grackle is a physical specific energy.

Light speed.
  ``ISRF_packet_light_speed_fraction_of_c`` sets a physical speed, so a share covers a comoving distance of :math:`v\,\Delta t / a` in a physical time :math:`\Delta t`. With ``fixed`` the bound on the time step uses the physical smoothing length :math:`a h`.

Redshift of the light in transit.
  Light loses energy to the expansion on its way. A band edge is fixed in physical energy, so the energy of a band falls, and the part of the Lyman-Werner light that drops below 11.2 eV enters the photoelectric band. The photon numbers change at their own rates. The loss of a hop grows with the physical path times :math:`H / c`: it is of the order of :math:`10^{-8}` for a 7 pc hop at :math:`z = 0` and of :math:`10^{-7}` at :math:`z = 9`.

Booking.
  The loss is taken once, when a share is booked at its receiver, and before the dust absorbs. The energy removed is booked as redshifted, not as absorbed, and the light that remains is what the dust sees. With the debugging option, the redshifted energy is a snapshot field of the gas (``PacketRedshiftedEnergies``). The photoelectric band gains light from the Lyman-Werner band, so its redshifted energy can be negative. Light that escapes (no downstream neighbour) is booked as escaped, as in a static box.

Accuracy.
  The shares are stored in single precision, so below :math:`z \approx 4` the loss of one hop is smaller than the rounding of a stored share, and the stored light carries it only on average. The ledgers are in double precision and book the loss exactly. A check of the energy balance must allow for this rounding.

Which mode?
~~~~~~~~~~~

Use the default (no transport) when only the gas inside the source kernels matters, for instance for photoelectric heating in a well-resolved interstellar medium. Use ``ISRF_propagation: 1`` when the light must reach gas beyond the kernels of its stars, be shadowed by the gas in between, or arrive after a delay. The transport costs memory in every gas particle (see below). Measured on one fixture (32768 gas particles in a periodic box, one star, all gas in one time bin, 12 slots, one thread, Grackle cooling on), a step with the transport costs 1.21 times a step without it during the first steps (steps 2 to 8, 904 against 749 ms). The cost then grows as the light fills the box and is 3.92 times the cost without transport from step 40 on (2946 against 752 ms). The cost follows the number of pairs of a particle and a non-empty slot, so a run with many stars reaches this plateau at once, and a fixture with more expensive cooling shows a smaller ratio. Runs with MPI, with several time bins or with many stars were not measured. With ``ISRF_packet_light_speed: fixed`` the time step of the gas is also bounded, which can make the run slower.

.. _gear_isrf_glossary:

Glossary
~~~~~~~~

Packet
  Energy of one band that travels from a source to a gas particle with one direction.
Emission cone
  Solid angle into which a star sends one share of its luminosity.
Reception slot
  One of a fixed number of cells of equal solid angle around a gas particle, fixed in a rotated frame. All light that arrives from one cell is merged into one packet.
Arrival bin, far cell
  Time cells of each slot. A near arrival bin is one time step wide. A far cell holds light that arrives later than the near bins.
Inbox
  Light received by a particle since its last step.
Reservoir
  Light that waits beyond the arrival bins.
Downstream neighbour
  For one reception slot of a gas particle: a gas neighbour inside the smoothing kernel of the particle (distance below :math:`\gamma_K h`) that lies ahead of the direction of the slot, with the cosine of the angle to that direction above :math:`10^{-3}`. Only these neighbours receive light from that slot. A slot without one lets its light escape.
Kernel mean
  The mean of the field of the neighbours of a particle, weighted by mass times the smoothing kernel :math:`W(r/h)`, with the particle itself included. This is the field that Grackle reads.
Ledger
  Snapshot fields that close the energy budget of a run (see :ref:`the packet fields <gear_output_isrf_packet>`).

Compile-time options
--------------------

``--with-isrf-packet-slots=<N>``
  Number of reception slots of each gas particle: 12 (default), 8 or 16. More slots keep the directions of incoming light apart better and cost memory in every gas particle. See :ref:`gear_radiation` for the shapes of the partitions.

``RADIATION_PACKET_BINS`` (default 4) and ``RADIATION_PACKET_FAR_CELLS`` (fixed at 2)
  Number of near arrival bins per slot, one time step wide, and of far arrival cells per slot, for light that arrives later. Set the first with ``CFLAGS="-DRADIATION_PACKET_BINS=<n>"`` at configure time. The defaults are enough for ``cap`` and ``fixed``.

``--with-chealpix=<path>``
  Path of the HEALPix C library (https://healpix.sourceforge.io), needed for the emission cones of the stars. HEALPix can need CFITSIO. See :doc:`the compilation instructions </GettingStarted/compiling_code>`.

Memory
~~~~~~

The transport adds 200 bytes per reception slot to every gas particle at the default number of bins. With ``--enable-isrf-debugging-checks`` off (the default), a gas particle is 2976 bytes with 12 slots and 2176 bytes with 8 slots (``sizeof(struct part)`` of a GEAR build with ``grackle_0``); with the option on, 3072 and 2272 bytes. The 16-slot build adds 800 bytes to the 12-slot size. For example, :math:`10^7` gas particles need 30 GB for the particles alone with 12 slots and 22 GB with 8 slots, without the option. A snapshot can be large too, see :ref:`gear_output_isrf`.

Restart and MPI
---------------

A restart file carries a layout marker. SWIFT refuses a restart file written by another code version, another number of slots or another ``--with-subgrid-radiation`` selection. Rerun from the initial conditions instead of resuming across versions.

The transport works with MPI: a star of a foreign rank hands its light to the local receivers. MPI support of the radiation model is still being tested, so check a multi-rank run against a single-rank one.

Limits
------

The light moves at most at the chosen speed in the frame of the gas, and no faster than :math:`c`. The transport is a hop between neighbouring gas particles, not a ray trace. Light that crosses a region of few particles spreads laterally by the transmission cone, and each reception slot merges all the light that arrives from one cell of directions into one packet with one direction.

Grackle coupling
----------------

Three ``GrackleCooling`` parameters control how the field is turned into heating and dissociation rates. Their defaults are chosen for a run without the ISRF, so a new ISRF run should look at all three.

``GrackleCooling:H2_self_shielding`` (default ``0``)
  How :math:`\mathrm{H}_2` shields itself from the Lyman-Werner field. ``0`` means no shielding, ``2`` sets the shielding column from the kernel support radius (see ``H2_self_shielding_path``), and ``3`` uses Grackle's own local Jeans length. Mode ``1`` is rejected: its length is read from neighbouring Cartesian grid cells, and SWIFT calls Grackle one particle at a time.

  **The default of 0 is the wrong choice for an ISRF run that tracks** :math:`\mathbf{H}_2`. Unshielded, the local Lyman-Werner rate dissociates molecular gas that a real cloud would keep. Set ``3`` unless you have a reason to prefer ``2``; the shipped ``ISRFH2Photodissociation`` example sets ``3``. SWIFT emits a start-up warning if the module is on, :math:`\mathrm{H}_2` is tracked and this is still ``0``.

``GrackleCooling:H2_self_shielding_path`` (default ``kernel_radius``)
  Used by mode ``2`` only. ``kernel_radius`` sets the shielding column path to one kernel support radius, ``kernel_diameter`` to twice the kernel support radius. Any other value stops SWIFT at start-up.

``GrackleCooling:photoelectric_heating_efficiency`` (default ``constant``)
  Which photoelectric efficiency Grackle applies to the Habing field (the PE+Lyman-Werner sum described above): ``constant`` (a fixed efficiency of 0.05 with prefactor :math:`10^{-24}`, Wolfire et al. 1995 Eq. 1), ``wolfire1995`` (the electron-density-dependent efficiency of the same paper, Eq. 2), or ``density_dependent`` (:math:`\min(0.07, 0.0149\,n_\mathrm{H}^{0.235})` with prefactor :math:`1.3 \times 10^{-24}`, Smith 2026 Eqs. A1-A2; requires a Grackle build that provides it; SWIFT stops at start-up if the module is on and the Grackle it was built against does not). None of the shipped ISRF examples override this: they run at the default, ``constant``.

Two further parameters are worth checking:

- ``GrackleCooling:local_dust_to_gas_ratio`` sets the dust content that the photoelectric rate is proportional to. A negative value uses Grackle's own default.
- ``GrackleCooling:RT_H2_dissociation_rate_cgs`` applies a constant :math:`\mathrm{H}_2` dissociation rate to all gas, on top of the per-particle rate this module computes. Set it to ``0`` in an ISRF run. SWIFT warns if it is nonzero.

Receiver-side extinction
------------------------

``GEARFeedback:ISRF_extinction_path`` (default ``pair_separation``) selects the mechanism that sets the path length :math:`l` of the dust column each gas particle shields itself with. The column is :math:`\Sigma = \rho_j l`, with the receiver's own density.

``constant_kernel_path``
  :math:`l = R\,\gamma_K h_j`, with :math:`R` set by ``ISRF_extinction_path_in_kernel_radii`` (default ``1.0``). :math:`R = 1` is one kernel support radius, the largest path the geometry admits, since the illuminating star sits inside the receiver's own kernel.

``pair_separation``
  :math:`l = r`, the star-to-particle separation of the pair being injected. It is the only mechanism that varies the attenuation across the kernel instead of applying one flat factor, and it is the default. Its kernel-weighted mean is :math:`5/12\,\gamma_K h` in the uniform optically thin limit, the value :math:`R = 5/12` of ``constant_kernel_path`` reproduces there.

``temperature_capped_jeans``
  :math:`l = \min(\lambda_J(\min(T, T_\mathrm{cap})), \gamma_K h_j)`, with :math:`T_\mathrm{cap}` set by ``ISRF_extinction_jeans_temperature_cap_K`` (default ``40`` K). Safranek-Shrader et al. (2017), MNRAS 465, 885, Section 3.4, rank the temperature-capped Jeans length best of five local column estimators. The cap at one support radius is not optional: uncapped, the Jeans length is orders of magnitude too opaque in diffuse gas, since the extinction acts only inside the illuminating star's own kernel.

An unrecognised value is a fatal error, not a fallback. This parameter is independent of ``GrackleCooling:H2_self_shielding_path``, since the two shield different processes.

Complete parameter list
-----------------------

The ISRF section of the ``GEARFeedback`` block, with every parameter at its default:

.. code:: YAML

   GEARFeedback:
     with_interstellar_radiation_field: 0                    # Master switch of the ISRF module
     ISRF_propagation: 0                                     # Transport the light with the packet hop
     ISRF_extinction_path: pair_separation                   # Receiver-side dust column path mechanism
     ISRF_extinction_path_in_kernel_radii: 1.0               # Path R in kernel support radii, constant_kernel_path only
     ISRF_extinction_jeans_temperature_cap_K: 40             # Jeans-length temperature cap, temperature_capped_jeans only
     ISRF_packet_light_speed_time_step_margin: 0.5           # Margin of the time step bound, ISRF_packet_light_speed fixed only
     ISRF_packet_light_speed: cap                            # Speed limit of the light: cap, fixed or unlimited
     ISRF_packet_light_speed_fraction_of_c: 0.99             # Limit speed as a fraction of c
     ISRF_packet_cone_power: 4                               # Power of the transmission cone
     ISRF_packet_read_mean_time_Myr: 0                       # Running mean time of the field read, 0 is off
     ISRF_packet_octant_frame_angles_deg: [24, 36, 66]       # Rotation of the frame of the reception slots
     ISRF_packet_octant_frame_seed_offset: 0                 # Seed of an extra random rotation of the frame

and the recommended ``GrackleCooling`` entries for an ISRF run tracking :math:`\mathrm{H}_2` (these are not Grackle's own defaults; see :ref:`gear_grackle_cooling` for the full block and its actual defaults):

.. code:: YAML

   GrackleCooling:
     H2_self_shielding: 3                         # Use 3 (local Jeans length) or 2 in an ISRF run tracking H2
     H2_self_shielding_path: kernel_radius        # Mode 2 only: kernel_radius or kernel_diameter
     photoelectric_heating_efficiency: constant   # constant, wolfire1995 or density_dependent
     local_dust_to_gas_ratio: -1                  # -1 uses Grackle's own default
     RT_H2_dissociation_rate_cgs: 0               # Keep at 0 in an ISRF run

Initial conditions
------------------

Three optional gas fields let an initial-conditions file seed the radiation field directly, for a test setup that starts with a field already in place rather than with a star that produces one. All three are singular, following SWIFT's convention that distinguishes input from output fields, and all three are for test and diagnostic use, not a normal production IC. With ``ISRF_propagation: 1`` the field is rewritten from the packet at every step.

.. list-table::
   :header-rows: 1
   :widths: 30 45 25

   * - Name
     - Description
     - Units
   * - ``PESpecificEnergy``
     - Initial specific PE-band energy
     - [U_L^2 U_T^{-2}]
   * - ``LWSpecificEnergy``
     - Initial specific Lyman-Werner-band energy
     - [U_L^2 U_T^{-2}]
   * - ``LWPhotonSpecificEnergy``
     - Initial Lyman-Werner-band photon-number moment, energy-equivalent at a fixed reference photon energy, not a photon count
     - [U_L^2 U_T^{-2}]

An initial-conditions file without any of them is unaffected: every field starts at zero. ``LWSpecificEnergy`` and ``LWPhotonSpecificEnergy`` are coupled at first init, per particle: on any particle where one of the two is zero and the other is not, SWIFT sets the zero one equal to the other, since a Lyman-Werner field with no accompanying photon moment is by definition at the reference photon energy. A particle with a nonzero ``LWSpecificEnergy`` therefore always carries a photon moment: to start with no photon moment, set ``LWSpecificEnergy`` to zero as well.

Snapshot outputs
------------------

The star and gas fields this channel writes are documented on the :ref:`gear_output_isrf` section of the :ref:`gear_output_fields` page. Only two gas fields are needed for an analysis: ``PESpecificEnergies`` and ``LWSpecificEnergies``. ``LWPhotonSpecificEnergies`` is the Lyman-Werner photon moment; no heating or dissociation rate reads it. The ``Packet*`` fields are energy ledgers and diagnostics, written only by a build with ``--enable-isrf-debugging-checks``. Switch them off with the output selection of :ref:`Output_selection_label` to save space.
