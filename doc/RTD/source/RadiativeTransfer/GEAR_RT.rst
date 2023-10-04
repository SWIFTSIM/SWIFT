.. GEAR Radiative Transfer
    Mladen Ivkovic 05.2021

.. _rt_GEAR:
   
GEAR RT
-------

.. warning::
    The radiative transfer schemes are still in development and are not useable
    at this moment. This page is currently a placeholder to document new
    features and requirements as the code grows.


Compiling for GEAR RT
~~~~~~~~~~~~~~~~~~~~~

-   To compile swift to be able to run with GEAR RT, you need to configure with
    ``--with-rt=GEAR_N`` where ``N`` is the integer number of photon groups that 
    you intend to use in your simulation.

-   You need to choose a Riemann solver for the RT equations. You can choose
    between the ``GLF`` and ``HLL`` solver. For the time being, I recommend 
    sticking to the ``GLF`` solver as the ``HLL`` solver is more expensive,
    but seemingly offers no advantage, although this remains to be confirmed
    in further testing.

-   GEAR RT is only compatible with the Meshless Finite Volume scheme. You'll
    need to compile using ``--with-hydro=gizmo-mfv``, which will also require
    you to select a hydro Riemann solver, e.g ``--with-riemann-solver=hllc``.

-   The thermochemistry requires the `grackle <https://github.com/grackle-project/grackle>`_ 
    library version above 3.2.1. [#f4]_ Grackle is a chemistry and cooling library presented in 
    `B. Smith et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`_.
    Please note that the current implementation is not (yet) as
    advanced as the :ref:`GEAR subgrid model grackle cooling <gear_grackle_cooling>`, 
    and the parameters listed as available there are not applicable for the 
    grackle cooling in combination with GEAR RT. You can however follow the Grackle 
    installation instructions documented there.

.. warning::
    (State 2023) Grackle is experiencing current development, and the API is subject
    to changes in the future. For convenience, a frozen version is hosted as a fork
    on github here: https://github.com/mladenivkovic/grackle-swift .
    The version available there will be tried and tested and ensured to work with
    GEAR-RT. 

    Additionally, that repository hosts files necessary to install that specific 
    version of grackle with spack.




Compulsory Runtime Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You need to provide the following runtime parameters in the yaml file (which will 
be further explained below):

.. code:: yaml

   GEARRT:
       photon_groups_Hz: [3.288e15, 5.945e15, 13.157e15]  # Photon frequency group bin edges in Hz
       stellar_spectrum_type: 0                           # Which radiation spectrum to use. 
                                                          #   0: constant. 
                                                          #   1: blackbody spectrum.
       stellar_luminosity_model: const                    # Which luminosity model to use.
       const_stellar_luminosities_LSol: [1., 1., 1.]      # stellar emission rates for each photon 
                                                          #   frequency bin in units of solar luminosity
                                                          #   for the 'const' luminosity model
       f_reduce_c: 1e-3                                   # reduce the speed of light by this factor
       CFL_condition: 0.9                                 # CFL condition for time integration
       hydrogen_mass_fraction:  0.76                      # total hydrogen (H + H+) mass fraction in the 
                                                          # metal-free portion of the gas


   TimeIntegration:
       max_nr_rt_subcycles: 128         # maximal number of RT subcycles per hydro step


The ``photon_groups_Hz`` need to be ``N`` frequency edges (floats) to separate 
the spectrum into ``N`` groups, where ``N`` is the same number you configured
with using ``--with_rt=GEAR_N``. The edges are **lower** edges of the bins, and
need to be sorted in increasing order. The final upper edge is defined in a 
different manner, and depends on the stellar spectrum type you assume (see below
for more details).

To specify the radiation emitted by stars, there are two main parameters:
``stellar_luminosity_model`` defines which model to use to obtain star 
luminosities, while ``stellar_spectrum_type`` determines the spectrum of the
radiation.
At the moment, the only way to define star emission rates is to use constant
stellar luminosities by setting ``stellar_luminosity_model: const``. [#f3]_
The constant star emission rates need to be provided in the parameter file and
to be defined for each photon frequency group individually using the 
``const_stellar_luminosities_LSol`` parameter. The luminosities are expected to
be in units of solar luminosities. Each star particle will then emit the given 
luminosities, independent of their other properties, e.g. the stellar age, 
metallicity, redshift, etc.

When solving the thermochemistry, we need to assume some form of stellar
spectrum so we may integrate over frequency bins to obtain average interaction
rates. The parameter ``stellar_spectrum_type`` is hence required, and allows you
to select between:

- constant spectrum (``stellar_spectrum_type: 0``)
    - Assume same energy density for any frequency.
    - This choice additionally requires you to provide a maximal frequency for
      the spectrum after which it'll be cut off via the 
      ``stellar_spectrum_const_max_frequency_Hz`` parameter

- blackbody spectrum (``stellar_spectrum_type: 1``)
    - Assume the spectrum is a blackbody spectrum
    - In this case, you need to provide also temperature of the blackbody via the 
      ``stellar_spectrum_blackbody_temperature_K`` parameter.
    - The assumed maximal considered frequency :math:`\nu_{max}` for this spectrum 
      is equal to 10 times :math:`\nu_{peak}`, the frequency at which the blackbody 
      spectrum has its maximum, i.e.

.. math::

     \nu_{peak} = 2.82144 \times k_{B} \times T / h_{Planck}

     \nu_{max} = 10 \times \nu_{peak}


.. warning::
   The ``stellar_spectrum_type`` parameter also determines the averaged photon 
   interaction cross sections, as they are being computed by integrating a 
   parametrization of the cross section multiplied by the assumed spectrum. See
   e.g. equations 9 - 11 in `Rosdahl et al. 2013. 
   <https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R/abstract>`_

Finally, you will also need to provide an upper threshold for the number of 
RT-subcycles w.r.t. a single hydro step via ``TimeIntegration:max_nr_rt_subcycles``.
For more details, refer to :ref:`the subcycling documentation <rt_subcycling>`.






Optional Runtime Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several optional parameters which can be set in the ``GEARRT`` block in the
yaml parameter file:

- ``f_limit_cooling_time``: (Default: ``0.6``) The cooling time computed by grackle 
  will be  multipled by this factor when estimating the particle's RT time step 
  size. If set to ``0.0``, the computation of cooling time is turned off.
- ``skip_thermochemistry``: (Default: ``0``) If set to ``1``, the entire 
  thermochemistry part of radiative transfer will be skipped. This is intended 
  only for debugging and testing the radiation transport, as it breaks the 
  purpose of RT.
- ``stars_max_timestep``: (Default: ``-1.``) Restrict the maximal timestep of stars 
  to this value (in internal units). Set to negative to turn off.
- ``grackle_verbose``: (Default: ``0``) set grackle to verbose.
- ``case_B_recombination``: (Default: ``1``) If ``1``, use case B recombination rates. 
  Otherwise, case A recombination rates are used.
- ``max_tchem_recursion``: (Default: ``0``) If set to positive nonzero value, the
  thermochemistry is computed using a "10% rule" *in addition to the ones grackle
  already uses*. If during a thermochemistry step the internal energy :math:``u`` of
  the gas changes by more than 10%, i.e. :math:`|u_{new}/u_{old} - 1| > 0.1`, then
  the thermochemistry step is repeated twice using half the time step size. This
  parameter sets the maximal recursion depth of halving the time step size and 
  repeating the entire thermochemistry step.

There are some further optional parameters related to setting up initial ion mass
fractions which are detailed in the 
:ref:`corresponding section of this documentation <rt_GEAR_set_ion_mass_fractions>`.










Choice of Internal Units
~~~~~~~~~~~~~~~~~~~~~~~~~~

The choice of internal units requires a bit of special attention. Part of the 
reason is that the exponents of the gas and radiation variables can quickly 
change by several dozens and cause overflows and other errors. Furthermore, the 
grackle library may have some other troubles with the units, e.g. when trying to
find a converging solution. [#f2]_

For this reason, I **strongly encourage** you to run the Internal Units check for 
GEAR-RT which you can find in the 
`swiftsim-rt-tools <https://github.com/SWIFTSIM/swiftsim-rt-tools/GEARRTUnitCheck>`_ 
repository under ``/GEARRTUnitsCheck``. The test should take no more than a 
minute to run, and requires only two yaml parameter files: the yaml parameter 
file that you intend to run your simulation with, and one that a provided script 
can extract automatically from the initial conditions hdf5 file. This test can 
save you a lot of headaches down the line.





Initial Conditions
~~~~~~~~~~~~~~~~~~

Setting Up Initial Conditions for RT
````````````````````````````````````

Optionally, you may want to provide initial conditions for the radiation field
and/or the mass fraction of the ionizing species.
To do so, you need to add the following datasets to the ``/PartType0`` particle
group:

.. code:: 

   PhotonEnergiesGroup1
   PhotonEnergiesGroup2 
   .
   .
   .
   PhotonEnergiesGroupN
   PhotonFluxesGroup1
   PhotonFluxesGroup2
   .
   .
   .
   PhotonFluxesGroupN
   MassFractionHI
   MassFractionHII
   MassFractionHeI
   MassFractionHeII
   MassFractionHeIII


-   The ``PhotonEnergies*`` datasets need to have dimension ``nparts``, while the
    ``PhotonFluxesGroup*`` datasets need to have dimension ``(nparts, 3)``, where
    ``nparts`` is the number of hydro particles. 
-   Note that the GEAR-RT scheme expects the ``PhotonEnergies*`` to be total 
    energies, not energy densities. 
-   If you are writing initial conditions where the fields have units [#f1]_, then 
    ``PhotonEnergies*`` are expected to have units of energy 
    :math:`[M L^2 T^{-2}]`), while the ``PhotonFluxes*`` fields should be in units 
    of energy times velocity (i.e. energy per unit time per unit area times volume, 
    :math:`[M L^3 T^{-3}]`).
-   The ``MassFraction*`` datasets need to have dimension ``nparts`` as well, and
    are all unitless.


Example using Python and ``swiftsimio``
````````````````````````````````````````

If you are using `swiftsimio <https://github.com/SWIFTSIM/swiftsimio>`_ to write
the initial condition files, then the easiest way of adding the RT initial
conditions is to first use the swiftsimio routines to write a file, then open it
up again and write the additional RT fields again using ``h5py`` routines.

Here is an example:

.. code:: python

    from swiftsimio import Writer
    import unyt
    import numpy as np
    import h5py

    # define unit system to use.
    unitsystem = unyt.unit_systems.cgs_unit_system

    # number of photon groups
    nPhotonGroups = 4

    # filename of ICs to be generated
    outputfilename = "my_rt_ICs.hdf5"

    # open a swiftsimio.Writer object
    w = Writer(...)

    # do your IC setup for gas, gravity etc now
    # ... 

    # write the IC file without doing anything RT related.
    w.write(outputfilename)

    # Now open file back up again and add RT data.
    F = h5py.File(outputfilename, "r+")
    header = F["Header"]
    nparts = header.attrs["NumPart_ThisFile"][0]
    parts = F["/PartType0"]

    # Create initial photon energies and fluxes. You can leave them unitless, 
    # the units have already been written down with w.write(). In this case, 
    # it's in cgs.
    for grp in range(nPhotonGroups):
        dsetname = "PhotonEnergiesGroup{0:d}".format(grp + 1)
        energydata = np.ones((nparts), dtype=np.float32) * some_value_you_want
        parts.create_dataset(dsetname, data=energydata)

        dsetname = "PhotonFluxesGroup{0:d}".format(grp + 1)
        fluxdata = np.zeros((nparts, 3), dtype=np.float32) * some_value_you_want
        parts.create_dataset(dsetname, data=fluxdata)

    # Create initial ionization species mass fractions.     
    HIdata = np.ones((nparts), dtype=np.float32) * 0.4
    parts.create_dataset("MassFractionHI", data=HIdata)
    HIIdata = np.ones((nparts), dtype=np.float32) * 0.1
    parts.create_dataset("MassFractionHII", data=HIIdata)
    HeIdata = np.ones((nparts), dtype=np.float32) * 0.3
    parts.create_dataset("MassFractionHeI", data=HeIdata)
    HeIIdata = np.ones((nparts), dtype=np.float32) * 0.15
    parts.create_dataset("MassFractionHeII", data=HeIIdata)
    HeIIIdata = np.ones((nparts), dtype=np.float32) * 0.05
    parts.create_dataset("MassFractionHeIII", data=HeIIIdata)

    # close up, and we're done!
    F.close()



.. _rt_GEAR_set_ion_mass_fractions:

Generate Ionization Mass Fractions Using SWIFT
``````````````````````````````````````````````

.. warning:: Using SWIFT to generate initial ionization mass fractions will
   overwrite the mass fractions that have been read in from the initial 
   conditions.

Optionally, you can use SWIFT to generate the initial mass fractions of the
ionizing species. To set the initial mass fractions of all particles to the same
value, use the following parameters in the yaml parameter file:

.. code:: yaml

    set_initial_ionization_mass_fractions: 1    # (Optional) manually overwrite initial mass fractions 
                                                # (using the values you set below)
    mass_fraction_HI: 0.76                      # set initial HI mass fractions to this value
    mass_fraction_HII: 0.                       # set initial HII mass fractions to this value
    mass_fraction_HeI: 0.24                     # set initial HeI mass fractions to this value
    mass_fraction_HeII: 0.                      # set initial HeII mass fractions to this value
    mass_fraction_HeIII: 0.                     # set initial HeIII mass fractions to this value

Alternatively, you can make SWIFT compute the initial ionization mass fractions
for you assuming ionization equilibrium, following `Katz, et al. 1996 
<ui.adsabs.harvard.edu/abs/1996ApJS..105...19K>`_ by setting

.. code:: yaml

    set_equilibrium_initial_ionization_mass_fractions: 1    # (Optional) set the initial ionization fractions 
                                                            # depending on gas temperature assuming ionization 
                                                            # equilibrium.
    hydrogen_mass_fraction:  0.76                           # total hydrogen (H + H+) mass fraction in the 
                                                            # metal-free portion of the gas

The ``hydrogen_mass_fraction`` (which is a compulsory argument in any case) will
determine the hydrogen and helium mass fractions, while SWIFT will determine the
equilibrium ionizations.

.. warning:: If you have somewhat sophisticated initial conditions (e.g. proper galaxies 
   etc) it is strongly recommended to set up the initial mass fractions to equilibrium
   values. Otherwise, the initial states can be too far off for GRACKLE's thermochemistry
   to handle it well, leading to all sorts of troubles, including crashes.




Accessing Output Data
~~~~~~~~~~~~~~~~~~~~~~

We recommend using `swiftsimio <https://github.com/SWIFTSIM/swiftsimio>`_ to 
access the RT related snapshot data. The compatibility is being maintained.
Here's an example how to access some specific quantities that you might find
useful:


.. code:: python

    #!/usr/bin/env python3

    import swiftsimio
    import unyt

    data = swiftsimio.load("output_0001.hdf5")
    meta = data.metadata



    # Accessing RT Related Metadata
    # ---------------------------------

    # get scheme name: "GEAR M1closure"
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    # number of photon groups used
    ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"])

    # get the reduced speed of light that was used. Will have unyts.
    reduced_speed_of_light = meta.reduced_lightspeed




    # Accessing Photon Data
    # ------------------------

    # accessing a photon group directly
    # NOTE: group names start with 1
    group_1_photon_energies = data.gas.photon_energies.group1
    group_1_photon_fluxes_x = data.gas.photon_fluxes.Group1X
    group_1_photon_fluxes_y = data.gas.photon_fluxes.Group1Y
    group_1_photon_fluxes_z = data.gas.photon_fluxes.Group1Z

    # want to stack all fluxes into 1 array?
    group1fluxes = swiftsimio.cosmo_array(
        unyt.uvstack(
            (group_1_photon_fluxes_x, group_1_photon_fluxes_y, group_1_photon_fluxes_z)
        ),
        group_1_photon_fluxes_x.units,
    ).T
    # group1fluxes.shape = (npart, 3)


    # Load all photon energies in a list
    photon_energies = [
        getattr(data.gas.photon_energies, "group" + str(g + 1)) for g in range(ngroups)
    ]



    # Accessing Ion Mass Fractions
    # -------------------------------
    fHI = data.gas.ion_mass_fractions.HI
    fHII = data.gas.ion_mass_fractions.HII
    fHeI = data.gas.ion_mass_fractions.HeI
    fHeII = data.gas.ion_mass_fractions.HeII
    fHeIII = data.gas.ion_mass_fractions.HeIII




.. rubric:: Footnotes

.. [#f1] To avoid possible confusions, here are some notes and equations
   regarding this choice of units.

   One of the RT equations solved by the GEAR RT is the zeroth moment of the
   equation of radiative transfer for each photon frequency group :math:`i` :

   :math:`\frac{\partial E_i}{\partial t} + \nabla \cdot \mathbf{F}_i = 0`

   where

   - :math:`E_i` : photon energy density; with :math:`[E_i] = erg / cm^3 = M L^{-1} T^{-2}`
   - :math:`F_i` : radiation flux (energy per unit time per unit surface); with :math:`[F_i] = erg / cm^2 / s = M T^{-3}` 

   and we neglect possible source and sink terms in this footnote.

   These dimensions are also used internally when solving the equations.
   For the initial conditions however, we require these quantities multiplied by
   the particle volume. The reason for this choice is so that the photon
   energies for each particle can be set by the users exactly, while the
   particle volume computation can be left to SWIFT to worry about internally.
   The addition of the particle volume term for the radiation flux was made so
   that the initial conditions are compatible with the SPHM1RT conventions, and
   both methods can run on the exact same ICs.


.. [#f2] For example, choosing cgs units as the internal units may lead to
   trouble with grackle. (Trouble like a gas at 10^6K without any heating
   sources heating up instead of cooling down.) The library is set up to work 
   with units geared towards cosmology. According to Britton Smith (private comm), 
   a decent rule of thumb is density_units ~ proton mass in g, time_units ~ 1 Myr 
   to 1 Gyr in s, length_units ~ 1 kpc to 1 Mpc in cm. This should keep you in a 
   relatively safe range.
   This is the state of things at 08.2022, with grackle being at version 3.2 (commit
   ``a089c837b8649c97b53ed3c51c84b1decf5073d8``)
    
.. [#f3] Technically there is also the model used for "Test 4" from the 
   `I. Iliev et al. 2006 <https://ui.adsabs.harvard.edu/abs/2006MNRAS.369.1625I>`_ 
   paper, but that is very specialized and shouldn't have much use in real 
   applications.

.. [#f4] Grackle version 3.2.1 still contained a bug related to the use of their
   "threadsafe functions" that could lead to disastrous outcomes. That was fixed
   in commit `a59489f`. So versions after 3.2.1 should work as expected.
   To be safe, we recommend you use the forked grackle repository specifically
   intended to freeze a stable version for use with SWIFT. You can find that fork
   on github: https://github.com/mladenivkovic/grackle-swift .
