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



Planetary EoS
-------------

See :ref:`planetary_eos`.



How to Implement a New Equation of State
----------------------------------------

See :ref:`new_option` for a full list of required changes.

You will need to provide an ``equation_of_state.h`` file containing: the
definition of ``eos_parameters``, IO functions and transformations between the
different variables: \\(u(\\rho, A)\\), \\(u(\\rho, P)\\), \\(P(\\rho,A)\\),
\\(P(\\rho, u)\\), \\(A(\\rho, P)\\), \\(A(\\rho, u)\\), \\(c_s(\\rho, A)\\),
\\(c_s(\\rho, u)\\) and \\(c_s(\\rho, P)\\). See other equation of state files
to have implementation details.
