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





Runtime Parameters
~~~~~~~~~~~~~~~~~~

You need to provide the following runtime parameters in the yaml file:

.. code:: yaml

   GEARRT:
       photon_groups_Hz: [3.288e15, 5.945e15, 13.157e15]  # Photon frequency group bin edges in Hz
       use_const_emission_rates: 1 
       star_emission_rates_LSol: [1., 1., 1., 1.]         # stellar emission rates for each photon 
                                                          # frequency bin in units of solar luminosity
       f_reduce_c: 1e-3                                   # reduce the speed of light by this factor
       CFL_condition: 0.99                                # CFL condition for time integration


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
intended to be optional in the future, for now it needs to be set to 1.



Initial Conditions
~~~~~~~~~~~~~~~~~~

Optionally, you may want to provide initial conditions for the radiation field.
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


The ``PhotonEnergiesX`` datasets need to have dimension ``nparts``, while the
``PhotonFluxesGroupX`` datasets need to have dimension ``(nparts, 3)``, where
``nparts`` is the number of hydro particles.

