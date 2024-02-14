.. AGN spin and jet model
   Filip Husko, 1 April 2022

.. AGN_spin_jet:

AGN jets and black hole spin in hydrodynamical simulations
==========================================================

Model summary
-------------


The main feature of this model, in terms of the effects on galaxy populations, is the addition of an AGN jet 
mode of feedback. In order to launch realistic jets, black hole spin is tracked and evolved for all BH 
particles in the simulation.

Jet powers, in addition to depending on spin, also depend on which accretion state the black hole is in. We 
include three accretion states: the thick, thin and slim disk. The thick disk appears at low accretion rates, 
has very strong jets and is inefficient at spinning up the black hole. The thin disk, appearing at 
intermediate accretion rates, has weak jets, strong radiation and efficiently spins up the black hole. The 
slim disk, corresponding to super-Eddington accretion, has features of both, and has both strong radiation and 
jets. Slim disks can be turned off in the model, but the thick and thin disks are intimitely tied to their 
feedback modes (jets and radiation, respectively).

In ``theory.rst`` we outline all of the theory which is implemented as part of the model. This includes when 
the black holes transition from one state to another, the strength of feedback in each state, how spin is 
evolved in terms of magnitude and direction, etc. In ``numerics.rst`` we discuss how jet launching is 
implemented, and additional black hole time steps introduced into the code. In ``params.rst`` we list and 
discuss all parameters used by the model. In ``output.rst`` we list additional arrays output for the BHs and 
tracers. Below we outline how to configure and run the model.

Compiling and running the model
-------------------------------

The model can be run with either the EAGLE or COLIBRE models. You can configure the model with ``--with-black-holes=SPIN_JET`` in combination with other configure options, or you can configure the full EAGLE or COLIBRE models with the new spin/jet physics as ``--with-subgrid=SPIN_JET_EAGLE`` and ``--with-subgrid=SPIN_JET_COLIBRE``, respectively. The model will then run as long as ``--black-holes`` is among the runtime options.

For cosmological simulations you do not need to do anything special, but for isolated runs (or any runs with black holes in the initial conditions), the ICs must include two new fields for all black holes: a scalar field representing black hole spins called ``Spins`` and a vector field representing the directions of the spin called ``AngularMomentumDirections``. The former should be between 0 and 1, while the latter should be normalized to 1.

A full list of all relevant parameters of the model is in ``params.rst``. We also briefly describe the most important parameters which need to be set to run the model, as well as how to run it in different configurations.

.. toctree::

  theory
  numerics
  params
  output
  variable_heating_temperatures
