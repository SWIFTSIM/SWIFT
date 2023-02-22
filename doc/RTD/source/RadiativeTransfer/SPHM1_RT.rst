.. SPHM1RT Radiative Transfer
    Tsang Keung Chan 01.2022

.. _rt_SPHM1:
   
SPHM1 RT
--------

SPHM1RT is the first two-moment radiative transfer on smoothed particle hydrodynamics (`Chan et al. 2021
<https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.5784C/abstract>`_). It solves the radiation energy and flux equations with a modified Eddington tensor closure. It is adaptive, efficient, and easy to parallelize.

.. warning::
    The radiative transfer schemes are still in development and are not useable
    at this moment. This page is currently a placeholder to document new
    features and requirements as the code grows.


Compiling for SPHM1-RT
~~~~~~~~~~~~~~~~~~~~~~

-   To compile swift to be able to run with SPHM1-RT, you need to configure with
    ``--with-rt=SPHM1RT_N`` where ``N`` is the integer number of photon groups that 
    you intend to use in your simulation. The number of photon groups for SPHM1RT has
    to be four for now (since it is hard-coded in thermo-chemistry solver).

-   SPHM1-RT is compatible with any SPH scheme. You'll
    need to compile using ``--with-hydro=sphenix`` or other SPH schemes, e.g. we have tested gadget2, minimal, and sphenix.

-   SPHM1-RT solves non-equilibrium with the `SUNDIALS <https://computing.llnl.gov/projects/sundials>` library, 
    which is SUite of Nonlinear and DIfferential/ALgebraic Equation Solvers. The SUNDIALS version has to be  5 . 
    You'll need to compile using ``--with-sundials=$SUNDIALS_ROOT``    
    SUNDIALS_ROOT is the root directory that contains the lib and include directories, e.g. on cosma:
    SUNDIALS_ROOT=/cosma/local/sundials/5.1.0/
    (IMPORTANT: Need SUNDIALS version  = 5). 
    The instructions of installing Sundials can be found, e.g., 
    `here <https://sundials.readthedocs.io/en/latest/Install_link.html>` or in this `website 
    <https://richings.bitbucket.io/chimes/user_guide/GettingStarted/sundials.html>`.



Runtime Parameters
~~~~~~~~~~~~~~~~~~

You need to provide the following runtime parameters in the yaml file:

.. code:: yaml

    SPHM1RT:
        cred: 2.99792458e10                                 # value of reduced speed of light for the RT solver in code unit
        CFL_condition: 0.1                                  # CFL condition for RT, independent of hydro 
        chi:  [0, 0, 0]                                     # (Optional) initial opacity in code unit for all gas particles
        photon_groups_Hz: [3.288e15, 5.945e15, 13.157e15]   # Photon frequency group bin edges in Hz.
        use_const_emission_rates: 1                         # (Optional) use constant emission rates for stars as defined with star_emission_rates_erg_LSol parameter
        star_emission_rates: [1e-32, 1e-32, 1e-32]          # (Optional) constant star emission rates (internal unit: energy/time) for each photon frequency group to use if use_constant_emission_rates is set.
        stellar_spectrum_type: 0                            # Which radiation spectrum to use. 0: constant from 0 until some max frequency set by stellar_spectrum_const_max_frequency_Hz. 1: blackbody spectrum.
        stellar_spectrum_const_max_frequency_Hz: 1.e17      # (Conditional) if stellar_spectrum_type=0, use this maximal frequency for the constant spectrum. 
        stars_max_timestep: -1.                             # (Optional) restrict the maximal timestep of stars to this value (in internal units). Set to negative to turn off.
        reinject:               1                           # (Optional) gather energy around injection radius and re-inject the energy


The ``photon_groups_Hz`` need to be ``N - 1`` frequency edges (floats) to separate 
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
``star_emission_rates`` parameter.

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

If ``reinject`` is set to 1, then the star will correct the radiation energy and 
flux within the injection radius and re-inject radiation. Combined with larger 
smoothing length, this can increase the isotropy of the radiation and the central 
radiation distribution. 



Thermo-chemistry parameters for RT
``````````````````````````````````
.. code:: yaml

    relativeTolerance:          1e-3                    # (Optional) Relative tolerance for SPHM1RT thermo-chemistry intergration
    absoluteTolerance:          1e-10                   # (Optional) Absolute tolerance for SPHM1RT thermo-chemistry integration
    explicitTolerance:          0.1                     # (Optional) Tolerance below which we use the explicit solution in SPHM1RT thermo-chemistry
    ionizing_photon_energy_erg: [3.0208e-11, 5.61973e-11, 1.05154e-10]  # (Optional) ionizing photon energy in erg averaged over frequency bins
    skip_thermochemistry: 0                             # (Optional) skip the thermochemistry. This is intended only for debugging and testing the radiation transport, as it breaks the purpose of RT.
    coolingon:              1                           # (Optional) switch for cooling (and photoheating), but photo-ionization will be ongoing even if coolingon==0 
    useparams:              1                           # (Optional) switch to use thermo-chemistry parameters from the parameter file
    sigma_cross:            [8.13e-18, 1e-32, 1e-32] # (Conditional) (if useparams=1) The cross section of ionizing photons for hydrogen (cm^2)
    alphaB:                 2.59e-13                    # (Conditional) (if useparams=1) The case B recombination coefficient for hydrogen (cgs)
    beta:                   3.1e-16                   # (Conditional) (if useparams=1) The collisional ionization coefficient for hydrogen (cgs)


ionizing_photon_energy_erg is the photon energy averaged over frequency within a frequency bin, given the radiation spectrum. In the default case, 
the first value corresponds to the bin from HI ionizing frequency to HeI ionizing frequency.
The second value is from HeI ionizing frequency to HeII ionizing frequency.
The third value is above HeII ionizing frequency.
The default values are calculated with T=1e5 K blackbody spectrum.
We currently assume the spectral shape is unchanged and universal, i.e. T=1e5 K blackbody spectrum everywhere. 
ionizing_photon_energy_erg is used to convert photon energy density to photon number density, photo-heating, and photo-ionization.



sigma_cross is also cross-section averaged within a frequency bin. 

Currently, SPHM1RT uses CVODE in SUNDIALS to solve non-equilibrium hydrogen and helium thermochemistry in three frequency bins,
from HI-HeII, HeII-HeIII and HeIII-inf. The precise coefficients will be published in Chan et al. in prep.,
but they can be found in src/rt_cooling_rates.h

Note that the first parameter in the thermo-chemistry array 
corresponds to the second parameter in injection array. For example, if
star_emission_rates: [0.0, 1.0, 0.0, 0.0], 
the star emits in the HI-HeII frequency and interacts with the first bin (8.13e-18):
sigma_cross:            [8.13e-18, 0.0, 0.0]

relativeTolerance, absoluteTolerance, and explicitTolerance are tolerances used in the CVODE calculation. 
These tolerances can be relaxed to increase the calculation speed, which could sacrifice accuracy.

We can also turn off thermochemistry or cooling for testing purpose by skip_thermochemistry and coolingon.
For testing purpose, we can also overwrite the thermo-chemistry parameters by setting useparams to 1
Currently, useparams==1 only works for pure hydrogen gas.




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


The ``PhotonEnergies*`` datasets need to have dimension ``nparts``, while the
``PhotonFluxesGroup*`` datasets need to have dimension ``(nparts, 3)``, where
``nparts`` is the number of hydro particles. If you are writing initial
conditions where the fields have units, then ``PhotonEnergies*`` are expected to
have units of energy :math:`[M L^2 T^{-2}]`), while the ``PhotonFluxes*`` fields
should be in units of energy times speed, :math:`[M L^3
T^{-3}]`).


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

    # Create initial element mass fractions.
    # Can be overwritten in parameter file if init_mass_fraction_metal is not -1.f (or set)
    # the element order: [Hydrogen, Helium]
    mfracH = np.ones(numPart)
    mfracHe = np.ones(numPart) * 0.0
    EMFdata = np.stack((mfracH, mfracHe), axis=-1)
    parts.create_dataset("RtElementMassFractions", data=EMFdata)

    # Create initial species abundances.
    # abundance is in n_X/n_H unit.
    # Can be overwritten in parameter file if useabundances = 1
    # the abundance order: [e, HI, HII, HeI, HeII, HeIII]
    Ae = np.ones(numPart) * 0.0   
    AHI = np.ones(numPart) * 1.0  
    AHII = np.ones(numPart) * 0.0 
    AHeI = np.ones(numPart) * 0.0 
    AHeII = np.ones(numPart) * 0.0 
    AHeIII = np.ones(numPart) * 0.0 
    SAdata = np.stack((Ae, AHI, AHII, AHeI, AHeII, AHeIII), axis=-1)    
    parts.create_dataset("RtSpeciesAbundances", data=SAdata)

    # close up, and we're done!
    F.close()



Generate Ionization Mass Fractions Using SWIFT
``````````````````````````````````````````````

.. warning:: Using SWIFT to generate initial ionization mass fractions will
   overwrite the mass fractions that have been read in from the initial 
   conditions.

Optionally, you can use SWIFT to generate the initial mass fractions of the
elements. To set the initial mass fractions of all particles to the same
value, use the following parameters in the yaml parameter file:

.. code:: yaml

  init_mass_fraction_metal:     0.                    # (Optional) Inital mass fraction of particle mass in *all* metals (if it is set or not equal to -1.F, the initial fraction will be over-written.)
  init_mass_fraction_Hydrogen:  1.0                   # (Conditional) (if init_mass_fraction_metal != -1.0f) Inital mass fraction of particle mass in Hydrogen
  init_mass_fraction_Helium:    0.0                   # (Conditional) (if init_mass_fraction_metal != -1.0f) Inital mass fraction of particle mass in Helium

To set the species abundances of all particles to the same
value, use the following parameters in the yaml parameter file:

.. code:: yaml

  useabundances:              1                       # (Optional) use the species abundances below, instead of reading from initial condition
  init_species_abundance_e:        1e-5               # (Conditional) (if useabundances==1) free electron abundances (in unit hydrogen number density:nH)
  init_species_abundance_HI:       0.99999            # (Conditional) (if useabundances==1) HI abundances (in unit hydrogen number density:nH)
  init_species_abundance_HII:      1e-5               # (Conditional) (if useabundances==1) HII abundances (in unit hydrogen number density:nH)
  init_species_abundance_HeI:      0.0                # (Conditional) (if useabundances==1) HeI abundances (in unit hydrogen number density:nH)
  init_species_abundance_HeII:     0.0                # (Conditional) (if useabundances==1) HeII abundances (in unit hydrogen number density:nH)
  init_species_abundance_HeIII:    0.0                # (Conditional) (if useabundances==1) HeIII abundances (in unit hydrogen number density:nH)


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

    # get scheme name: "SPH M1closure"
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


    # Accessing Element mass fraction
    fH = data.gas.rt_element_mass_fractions.hydrogen
    fHe = data.gas.rt_element_mass_fractions.helium

    # Accessing Species Abundances 
    # abundance is in n_X/n_H unit.
    # -------------------------------
    Ae = data.gas.rt_species_abundances.e
    AHI = data.gas.rt_species_abundances.HI
    AHII = data.gas.rt_species_abundances.HII
    AHeI = data.gas.rt_species_abundances.HeI
    AHeII = data.gas.rt_species_abundances.HeII
    AHeIII = data.gas.rt_species_abundances.HeIII
