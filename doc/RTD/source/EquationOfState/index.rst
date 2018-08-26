.. Equation of State
   Loic Hausammann, 6th April 2018

.. _equation_of_state:

Equation of State
=================

Currently (if the documentation was well updated), we have two different
equation of states implemented: ideal gas and isothermal.  They describe the
relations between our main thermodynamical variables: the internal energy
(\\(u\\)), the density (\\(\\rho\\)), the entropy (\\(A\\)) and the pressure
(\\(P\\)).

Equations
---------

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


How to Implement a New Equation of State
----------------------------------------

See :ref:`new_option` for a full list of required changes.

You will need to provide a ``equation_of_state.h`` file containing: the
definition of ``eos_parameters``, IO functions and transformations between the
different variables: \\(u(\\rho, A)\\), \\(u(\\rho, P)\\), \\(P(\\rho,A)\\),
\\(P(\\rho, u)\\), \\(A(\\rho, P)\\), \\(A(\\rho, u)\\), \\(c_s(\\rho, A)\\),
\\(c_s(\\rho, u)\\) and \\(c_s(\\rho, P)\\). See other equation of state files
to have implementation details.
