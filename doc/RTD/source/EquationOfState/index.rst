.. Equations of State
   Loic Hausammann, 6th April 2018
   Jacob Kegerreis, 13th March 2020

.. _equation_of_state:

Equations of State
==================

Currently, SWIFT offers three different gas equations of state (EoS)
implemented: ``ideal``, ``isothermal``, and ``barotropic``; as well as a variety
of EoS for "planetary" materials.  The EoS describe the relations between our
main thermodynamical variables: the internal energy per unit mass :math:`u`, the
mass density :math:`\rho`, the entropy :math:`A` and the pressure :math:`P`.
It is selected af configure time via the option ``--with-equation-of-state``.

Gas EoS
-------

We write the adiabatic index as :math:`\gamma` and :math:`c_s` denotes
the speed of sound. The adiabatic index can be changed at configure
time by choosing one of the allowed values of the option
``--with-adiabatic-index``. The default value is :math:`\gamma = 5/3`.

The tables below give the expression for the thermodynamic quantities
on each row entry as a function of the gas density and the
thermodynamical quantity given in the header of each column.

.. csv-table:: Ideal Gas
   :header: "Variable", "A", "u", "P"
	   
   "A", "", :math:`\left( \gamma - 1 \right) u \rho^{1-\gamma}`, :math:`P \rho^{-\gamma}`
   "u", :math:`A \frac{ \rho^{ \gamma - 1 } }{\gamma - 1 }`, "", :math:`\frac{1}{\gamma - 1} \frac{P}{\rho}`
   "P", :math:`A \rho^\gamma`, :math:`\left( \gamma - 1\right) u \rho`, ""
   :math:`c_s`, :math:`\sqrt{ \gamma \rho^{\gamma - 1} A}`, :math:`\sqrt{ u \gamma \left( \gamma - 1 \right) }`, :math:`\sqrt{ \frac{\gamma P}{\rho} }`


.. csv-table:: Isothermal Gas
   :header: "Variable", "-", "-", "-"

	    
   "A", "", :math:`\left( \gamma - 1 \right) u \rho^{1-\gamma}`, "" 
   "u", "", "const", ""
   "P", "", :math:`\left( \gamma - 1\right) u \rho`, ""
   :math:`c_s`, "", :math:`\sqrt{ u \gamma \left( \gamma - 1 \right) }`, ""

.. csv-table:: Barotropic Gas
   :header: "Variable", "-", "-", "-"

   "A", "", :math:`\rho^{1-\gamma} c_0^2 \sqrt{1 + \left( \frac{\rho}{\rho_c}  \right) }`, ""
   "u", "", :math:`\frac{1}(\gamma -1)c_0^2 \sqrt{1 + \left( \frac{\rho}{\rho_c}  \right) }`, ""
   "P", "", :math:`\rho c_0^2 \sqrt{1 + \left( \frac{\rho}{\rho_c}  \right) }`, ""
   :math:`c_s`, "", :math:`\sqrt{ c_0^2 \sqrt{1 + \left( \frac{\rho}{\rho_c}  \right) }}`, ""
   
Note that when running with an isothermal or barotropic equation of state, the
value of the tracked thermodynamic variable (e.g. the entropy in a
density-entropy scheme or the internal enegy in a density-energy SPH
formulation) written to the snapshots is meaningless. The pressure, however, is
always correct in all scheme.

For the isothermal equation of state, the internal energy is specified at
runtime via the parameter file. In the case of the barotropic gas, the vacuum
sound speed :math:`c_0` and core density :math:`\rho_c` are similarly specified.


Planetary EoS
-------------

See :ref:`planetary_eos`.



How to Implement a New Equation of State
----------------------------------------

See :ref:`new_option` for a full list of required changes.

You will need to provide an ``equation_of_state.h`` file containing: the
definition of ``eos_parameters``, IO functions and transformations between the
different variables: :math:`u(\rho, A)`, :math:`u(\rho, P)`, :math:`P(\rho,A)`,
:math:`P(\rho, u)`, :math:`A(\rho, P)`, :math:`A(\rho, u)`, :math:`c_s(\rho, A)`,
:math:`c_s(\rho, u)` and :math:`c_s(\rho, P)`. See other equation of state files
to have implementation details.
