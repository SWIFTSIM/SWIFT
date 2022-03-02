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
    but seemingly offers no advantage, although this remains to be comfirmed
    in further testing.

-   GEAR RT is only compatible with the Meshless Finite Volume scheme. You'll
    need to compile using ``--with-hydro=gizmo-mfv``, which will also require
    you to select a hydro Riemann solver, e.g ``--with-riemann-solver=hllc``.

-   The thermochemistry requires the `grackle <https://github.com/grackle-project/grackle>`_ 
    library. Grackle is a chemistry and cooling library presented in 
    `B. Smith et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`_.
    Please note that the current implementation is not (yet) as
    advanced as the :ref:`GEAR subgrid model grackle cooling <gear_grackle_cooling>`, 
    and the parameters listed as available there are not applicable for the 
    grackle cooling in combination with GEAR RT. 



Compulsory Runtime Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You need to provide the following runtime parameters in the yaml file:

.. code:: yaml

   GEARRT:
       photon_groups_Hz: [3.288e15, 5.945e15, 13.157e15]  # Photon frequency group bin edges in Hz
       use_const_emission_rates: 1 
       star_emission_rates_LSol: [1., 1., 1., 1.]         # stellar emission rates for each photon 
                                                          # frequency bin in units of solar luminosity
       f_reduce_c: 1e-3                                   # reduce the speed of light by this factor
       CFL_condition: 0.9                                 # CFL condition for time integration
       hydrogen_mass_fraction:  0.76                      # total hydrogen (H + H+) mass fraction in the 
                                                          # metal-free portion of the gas

       stellar_spectrum_type: 0                           # Which radiation spectrum to use. 0: constant. 1: blackbody spectrum.

The ``photon_groups`` need to be ``N - 1`` frequency edges (floats) to separate 
the spectrum into ``N`` groups. The outer limits of zero and infinity are 
assumed.

At the moment, the only way to define star emission rates is to use constant
star emission rates that need to be provided in the parameter file. The star 
emission rates need to be defined for each photon frequency group individually.
The first entry of the array is for the photon group with frequency 
``[0, <first entry of photon_groups_Hz>)``. Each star particle will then emit
the given energies, independent of their other properties.

Furthermore, even though the parameter ``use_const_emission_rates`` is 
intended to be optional in the future, **for now it needs to be set to 1.**, and
it requires you to manually set the stellar emission rates via the
``star_emission_rates_LSol`` parameter.

When solving the thermochemistry, we need to assume some form of stellar
spectrum so we may integrate over frequency bins to obtain average interaction
rates. The parameter ``stellar_spectrum_type`` is hence required, and allows you
to select between:

- constant spectrum (``stellar_spectrum_type: 0``)
    - This choice additionally requires you to provide a maximal frequency for
      the spectrum after which it'll be cut off via the 
      ``stellar_spectrum_const_max_frequency_Hz`` parameter

- blackbody spectrum (``stellar_spectrum_type: 1``)
    - In this case, you need to provide also temperature of the blackbody via the 
      ``stellar_spectrum_blackbody_temperature_K`` parameter.






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
-   If you are writing initial conditions where the fields have units, then 
    ``PhotonEnergies*`` are expected to have units of energy 
    :math:`[M L^2 T^{-2}]`), while the ``PhotonFluxes*`` fields should be in units 
    of energy flux (energy per unit time per unit area, :math:`[M T^{-3}]`). 
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
