.. RT Subcycling
    Mladen Ivkovic 07.2022

.. _rt_subcycling:
   
RT Subcycling
-------------

.. warning::
    The radiative transfer schemes are still in development and are not useable
    at this moment. This page is currently a placeholder to document new
    features and requirements as the code grows.


SWIFT allows to sub-cycle the solution of radiative transfer steps (both 
photon propagation and thermochemistry) with respect to the hydrodynamics
time steps. Basically you can tell SWIFT to run up to X radiative transfer
steps during a single hydrodynamics step for all particles in the simulation.
The aim is to not waste time doing unnecessary hydrodynamics updates, which
typically allow for much higher time steps compared to radiation due to the
propagation speed of the respective advected quantity.

You will need to provide an upper limit on how many RT subcycles per hydro
step you want to allow. That is governed by the

.. code:: yaml

   TimeIntegration:
       max_nr_rt_subcycles: 128         # maximal number of RT subcycles per hydro step

parameter, which is mandatory for any RT runs. To turn off subcycling and 
couple the radiative transfer and the hydrodynamics time steps one-to-one,
set this parameter to either 0 or 1.

Due to the discretization of individual particle time steps in time bins
with a factor of 2 difference in time step size from a lower to a higher
time bin, the ``max_nr_rt_subcycles`` parameter itself is required to be
a power of 2 as well.

Note that this parameter will set an upper limit to the number of subcycles
per hydro step. If the ratio of hydro-to-RT time step is greater than what
``max_nr_rt_subcycles`` allows for, then the hydro time step will be reduced
to fit the maximal threshold. If it is smaller, the particle will simply do 
fewer subcycles.

