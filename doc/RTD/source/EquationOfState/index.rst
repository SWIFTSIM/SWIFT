.. Equations of State
   Loic Hausammann, 6th April 2018
   Jacob Kegerreis, 13th March 2020

.. _equation_of_state:

Equations of State
==================

Currently, SWIFT offers two different gas equations of state (EoS)
implemented: ``ideal`` and ``isothermal``; as well as a variety of EoS for
"planetary" materials.  The EoS describe the relations between our
main thermodynamical variables: the internal energy per unit mass
(\\(u\\)), the mass density (\\(\\rho\\)), the entropy (\\(A\\)) and
the pressure (\\(P\\)).

Gas EoS
-------

We write the adiabatic index as \\(\\gamma \\) and \\( c_s \\) denotes
the speed of sound. The adiabatic index can be changed at configure
time by choosing one of the allowed values of the option
``--with-adiabatic-index``. The default value is \\(\\gamma = 5/3 \\).

The tables below give the expression for the thermodynamic quantities
on each row entry as a function of the gas density and the
thermodynamical quantity given in the header of each column.

.. csv-table:: Ideal Gas
   :header: "Variable", "A", "u", "P"
	   
   "A", "", "\\( \\left( \\gamma - 1 \\right) u \\rho^{1-\\gamma} \\)", "\\(P \\rho^{-\\gamma} \\)"
   "u", "\\( A \\frac{ \\rho^{ \\gamma - 1 } }{\\gamma - 1 } \\)", "", "\\(\\frac{1}{\\gamma - 1} \\frac{P}{\\rho}\\)"
   "P", "\\( A \\rho^\\gamma \\)", "\\( \\left( \\gamma - 1\\right) u \\rho \\)", ""
   "\\(c_s\\)", "\\(\\sqrt{ \\gamma \\rho^{\\gamma - 1} A}\\)", "\\(\\sqrt{ u \\gamma \\left( \\gamma - 1 \\right) } \\)", "\\(\\sqrt{ \\frac{\\gamma P}{\\rho} }\\)"


.. csv-table:: Isothermal Gas
   :header: "Variable", "A", "u", "P"

	    
   "A", "", "\\(\\left( \\gamma - 1 \\right) u \\rho^{1-\\gamma}\\)", "" 
   "u", "", "const", ""
   "P", "", "\\(\\left( \\gamma - 1\\right) u \\rho \\)", ""
   "\\( c_s\\)", "", "\\(\\sqrt{ u \\gamma \\left( \\gamma - 1 \\right) } \\)", ""

Note that when running with an isothermal equation of state, the value
of the tracked thermodynamic variable (e.g. the entropy in a
density-entropy scheme or the internal enegy in a density-energy SPH
formulation) written to the snapshots is meaningless. The pressure,
however, is always correct in all scheme.


.. _planetary_equation_of_state:

Planetary EoS
-------------
Configuring SWIFT with the ``--with-equation-of-state=planetary`` and 
``--with-hydro=planetary`` options enables the use of multiple EoS.
Every SPH particle then requires and carries the additional ``MaterialID`` flag 
from the initial conditions file. This flag indicates the particle's material 
and which EoS it should use. 

So far, we have implemented several Tillotson, SESAME, and Hubbard \& MacFarlane 
(1980) materials, with more on their way.
The material's ID is set by a base type ID (multiplied by 100), plus a minor 
type:

+ Tillotson (Melosh, 2007): Base type ``1``
    + Iron: ``100``
    + Granite: ``101``
    + Water: ``102``
+ Hubbard \& MacFarlane (1980): Base type ``2``
    + Hydrogen-helium atmosphere: ``200``
    + Ice H20-CH4-NH3 mix: ``201``
    + Rock SiO2-MgO-FeS-FeO mix: ``202``
+ SESAME (and similar): Base type ``3``
    + Iron (2140): ``300``
    + Basalt (7530): ``301``
    + Water (7154): ``302``
    + Senft \& Stewart (2008) water (in a SESAME-style table): ``303``

Unlike the EoS for an ideal or isothermal gas, these more complicated materials 
do not always include transformations between the internal energy, 
temperature, and entropy. At the moment, we have only implemented 
\\(P(\\rho, u)\\) and \\(c_s(\\rho, u)\\). 
This is sufficient for the simple :ref:`planetary_sph` hydrodynamics scheme, 
but makes these materials currently incompatible with entropy-based schemes.

The Tillotson sound speed was derived using 
\\(c_s^2 = \\left. \\dfrac{\\partial P}{\\partial \\rho} \\right|_S \\)
as described in `Kegerreis et al. (2019) <https://doi.org/10.1093/mnras/stz1606>`_.
The table files for the HM80 and SESAME-style EoS can be downloaded using 
the ``swiftsim/examples/EoSTables/get_eos_tables.sh`` script.

See :ref:`planetary` for other related information.


How to Implement a New Equation of State
----------------------------------------

See :ref:`new_option` for a full list of required changes.

You will need to provide an ``equation_of_state.h`` file containing: the
definition of ``eos_parameters``, IO functions and transformations between the
different variables: \\(u(\\rho, A)\\), \\(u(\\rho, P)\\), \\(P(\\rho,A)\\),
\\(P(\\rho, u)\\), \\(A(\\rho, P)\\), \\(A(\\rho, u)\\), \\(c_s(\\rho, A)\\),
\\(c_s(\\rho, u)\\) and \\(c_s(\\rho, P)\\). See other equation of state files
to have implementation details.
