.. Basic sub-grid model
   Matthieu Schaller, 20th December 2018


Basic model
===========


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
