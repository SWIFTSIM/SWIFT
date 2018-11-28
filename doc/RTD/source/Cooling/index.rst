.. Equation of State
   Loic Hausammann, 7th April 2018

.. _cooling:

Cooling
=======

Currently, we have 5 different cooling (EAGLE, Grackle, const-lambda, const-du
and none).  Three of them are easily solved analytically (const-lambda,
const-du and none) while the two last requires complex chemical networks.


Equations
---------

The first table compares the different analytical cooling while the next ones
are specific to a given cooling.  The quantities are the internal energy (\\( u
\\)), the density \\( rho \\), the element mass fraction (\\( X_i \\)), the
cooling function (\\(\\Lambda\\), the proton mass (\\( m_H \\)) and the time
step condition (\\( t\_\\text{step}\\)).  If not specified otherwise, all
cooling contains a temperature floor avoiding negative temperature.

.. csv-table:: Analytical Cooling
   :header: "Variable", "Const-Lambda", "Const-du", "None"

   "\\( \\frac{ \\mathrm{d}u }{ \\mathrm{d}t } \\)", "\\( -\\Lambda \\frac{\\rho^2 X_H^2}{\\rho m_H^2} \\)", "const", "0"
   "\\( \\Delta t\_\\text{max} \\)", "\\( t\_\\text{step} \\frac{u}{\\left|\\frac{ \\mathrm{d}u }{ \\mathrm{d}t }\\right|} \\)", "\\( t\_\\text{step} \\frac{u}{\\ \\left| \\frac{ \\mathrm{d}u }{ \\mathrm{d}t }\\right|} \\)", "None"


Grackle
~~~~~~~
   
Grackle is a chemistry and cooling library presented in `B. Smith et al. 2016 <https://arxiv.org/abs/1610.09591>`_ 
(do not forget to cite if used).  Four different modes are available:
equilibrium, 6 species network (H, H\\( ^+ \\), e\\( ^- \\), He, He\\( ^+ \\)
and He\\( ^{++} \\)), 9 species network (adds H\\(^-\\), H\\(_2\\) and
H\\(_2^+\\)) and 12 species (adds D, D\\(^+\\) and HD).  Following the same
order, the swift cooling options are ``grackle``, ``grackle1``, ``grackle2``
and ``grackle3`` (the numbers correspond to the value of
``primordial_chemistry`` in Grackle).  It also includes some self-shielding
methods and UV background.  In order to use the Grackle cooling, you will need
to provide an HDF5 table computed by Cloudy.

When starting a simulation without providing the different fractions, the code
supposes an equilibrium and computes the fractions automatically.

In order to compile SWIFT with Grackle, you need to provide the options ``with-grackle``
and ``with-chemistry``.

You will need a Grackle version later than 3.1. To compile it, run
the following commands from the root directory of Grackle:
``./configure; cd src/clib``.
Update the variables ``LOCAL_HDF5_INSTALL`` and ``MACH_INSTALL_PREFIX`` in
the file ``src/clib/Make.mach.linux-gnu``.
Finish with ``make machine-linux-gnu; make && make install``.
If you encounter any problem, you can look at the `Grackle documentation <https://grackle.readthedocs.io/en/latest/>`_

You can now provide the path given for ``MACH_INSTALL_PREFIX`` to ``with-grackle``.

Eagle
~~~~~

TODO

How to Implement a New Cooling
------------------------------

The developper should provide at least one function for:
 * writing the cooling name in HDF5
 * cooling a particle
 * the maximal time step possible
 * initializing a particle
 * computing the total energy radiated by a particle
 * initializing the cooling parameters
 * printing the cooling type

For implementation details, see ``src/cooling/none/cooling.h``

See :ref:`new_option` for the full list of changes required.
