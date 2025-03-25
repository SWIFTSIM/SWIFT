.. Basic sub-grid model
   Matthieu Schaller, 20th December 2018


Basic model (others)
====================

Sinks: Simple Bondi-Hoyle accretion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

The ``Basic`` sink model provides a foundation on which new sink implementations could be built. It includes a prescription for Bondi-Hoyle gas accretion, and a method for sink-sink mergers that is a slightly simplified version of the implementation used in GEAR.

No other physics is implemented for this model. Sinks cannot form - to use this model, sink particles must already be present in the initial conditions. They also cannot spawn stars from the gas they accrete.

Bondi-Hoyle accretion can be done in one of two ways:

 * Gas particles within the sink's kernel are stochastically swallowed entirely, with a probability set by the Bondi-Hoyle rate. Specifically, the probability is set by the current difference between the sink's subgrid mass (determined by the accretion rate) and its dynamical mass (which tracks the number of particles/sinks actually swallowed). This mode is equivalent to the EAGLE black hole accretion model.
 * Gas particles within the sink's kernel are "nibbled" down to some minimal mass, which can be specified by the user. This method is equivalent to the black hole accretion model of Bahe et al. 2022.

This model has only two parameters that must be specified in your parameter ``yml`` file:

 * ``BasicSink:use_nibbling``: determines whether accretion is done by "nibbling" or by swallowing outright.
 * ``BasicSink:min_gas_mass_for_nibbling_Msun``: if using "nibbling", the minimum mass to which gas particles can be nibbled. A good default is half the original particle mass.

For an even more bare-bones starting point, the ``Default`` sink model contains no physics at all, and is a totally blank canvas on which to build your sink model.


Cooling: Analytic models
~~~~~~~~~~~~~~~~~~~~~~~~

Currently, we have 3 different simple cooling models (const-lambda, const-du
and Compton). These are all based on analytic formulas and can be used
to quickly understand how the cooling interacts with the rest of the
code before moving to more complex models.

Equations
---------

The first table compares the different analytical cooling while the next ones
are specific to a given cooling.  The quantities are the internal energy (\\( u
\\)), the density \\( rho \\), the element mass fraction (\\( X_i \\)), the
cooling function (\\(\\Lambda\\), the proton mass (\\( m_H \\)) and the time
step condition (\\( t\_\\text{step}\\)).  If not specified otherwise, all
cooling contains a temperature floor avoiding negative temperature.

.. csv-table:: Analytical Cooling
   :header: "Variable", "Const-Lambda", "Const-du"

   "\\( \\frac{ \\mathrm{d}u }{ \\mathrm{d}t } \\)", "\\( -\\Lambda \\frac{\\rho^2 X_H^2}{\\rho m_H^2} \\)", "const"
   "\\( \\Delta t\_\\text{max} \\)", "\\( t\_\\text{step} \\frac{u}{\\left|\\frac{ \\mathrm{d}u }{ \\mathrm{d}t }\\right|} \\)", "\\( t\_\\text{step} \\frac{u}{\\ \\left| \\frac{ \\mathrm{d}u }{ \\mathrm{d}t }\\right|} \\)"

TODO: Add description of the parameters and units.

TODO: Add Compton cooling model

How to Implement a New Cooling
------------------------------

The developer should provide at least one function for:
 * writing the cooling name in HDF5
 * cooling a particle
 * the maximal time step possible
 * initializing a particle
 * computing the total energy radiated by a particle
 * initializing the cooling parameters
 * printing the cooling type

For implementation details, see ``src/cooling/none/cooling.h``

See :ref:`new_option` for the full list of changes required.
