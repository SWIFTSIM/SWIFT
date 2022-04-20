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
    you intend to use in your simulation.

-   SPHM1-RT is compatible with any SPH scheme. You'll
    need to compile using ``--with-hydro=sphenix`` or other SPH schemes, e.g. we have tested gadget2, minimal, and sphenix.




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
        star_emission_rates_LSol: [1e-32, 1e-32, 1e-32]     # (Optional) constant star emission rates for each photon frequency group to use if use_constant_emission_rates is set.
        stellar_spectrum_type: 0                            # Which radiation spectrum to use. 0: constant from 0 until some max frequency set by stellar_spectrum_const_max_frequency_Hz. 1: blackbody spectrum.
        stellar_spectrum_const_max_frequency_Hz: 1.e17      # (Conditional) if stellar_spectrum_type=0, use this maximal frequency for the constant spectrum. 
        stars_max_timestep: -1.                             # (Optional) restrict the maximal timestep of stars to this value (in internal units). Set to negative to turn off.


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

    # close up, and we're done!
    F.close()
