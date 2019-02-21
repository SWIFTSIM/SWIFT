.. Equations of State
   Loic Hausammann, 6th April 2018
   Jacob Kegerreis, 3rd February 2019

.. _equation_of_state:

Equations of State
==================

Currently (if the documentation was well updated), we have two different gas
equations of state (EoS) implemented: ideal and isothermal; as well as a variety  
of EoS for "planetary" materials. 
The EoS describe the relations between our main thermodynamical variables: 
the internal energy (\\(u\\)), the density (\\(\\rho\\)), the entropy (\\(A\\)) 
and the pressure (\\(P\\)).

Gas EoS
-------

In the following section, the variables not yet defined are: \\(\\gamma\\) for
the adiabatic index and \\( c_s \\) for the speed of sound.

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

+ Tillotson (Melosh, 2007): ``1``
    + Iron: ``100``
    + Granite: ``101``
    + Water: ``102``
+ Hubbard \& MacFarlane (1980): ``2``
    + Hydrogen-helium atmosphere: ``200``
    + Ice H20-CH4-NH3 mix: ``201``
    + Rock SiO2-MgO-FeS-FeO mix: ``202``
+ SESAME (and similar): ``3``
    + Iron (2140): ``300``
    + Basalt (7530): ``301``
    + Water (7154): ``302``
    + Senft \& Stewart (2008) water (in a SESAME-style table): ``303``

Unlike the EoS for an ideal or isothermal gas, these more complicated materials 
do not always include transformations between the internal energy, 
temperature, and entropy. At the moment, we have only implemented 
\\(P(\\rho, u)\\) and \\(c_s(\\rho, u)\\). 
This is sufficient for the simple :ref:`planetary_sph` hydrodynamics scheme, 
but makes these materials currently incompatible with other entropy-based 
schemes.

The Tillotson sound speed was derived using 
\\(c_s^2 = \\left. \\dfrac{\\partial P}{\\partial \\rho} \\right|_S \\)
as described in Kegerreis et al. (2019).
The table files for the HM80 and SESAME-style EoS can be downloaded using 
the ``examples/EoSTables/get_eos_tables.sh`` script.


How to Implement a New Equation of State
----------------------------------------

See :ref:`new_option` for a full list of required changes.

You will need to provide an ``equation_of_state.h`` file containing: the
definition of ``eos_parameters``, IO functions and transformations between the
different variables: \\(u(\\rho, A)\\), \\(u(\\rho, P)\\), \\(P(\\rho,A)\\),
\\(P(\\rho, u)\\), \\(A(\\rho, P)\\), \\(A(\\rho, u)\\), \\(c_s(\\rho, A)\\),
\\(c_s(\\rho, u)\\) and \\(c_s(\\rho, P)\\). See other equation of state files
to have implementation details.
