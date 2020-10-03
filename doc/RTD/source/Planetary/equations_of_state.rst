.. Planetary EoS
    Jacob Kegerreis, 13th March 2020

.. _planetary_eos:

Planetary Equations of State
============================
   
Configuring SWIFT with the ``--with-equation-of-state=planetary`` and 
``--with-hydro=planetary`` options enables the use of multiple 
equations of state (EoS).
Every SPH particle then requires and carries the additional ``MaterialID`` flag 
from the initial conditions file. This flag indicates the particle's material 
and which EoS it should use. 

It is important to check that the EoS you use are appropriate 
for the conditions in the simulation that you run.
Please follow the original sources of these EoS for more information and 
to check the regions of validity.

So far, we have implemented several Tillotson, ANEOS, SESAME, 
and Hubbard \& MacFarlane (1980) materials, with more on the way.
The material's ID is set by a base type ID (multiplied by 100), 
plus a minor type:

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
    + Senft \& Stewart (2008) water in a SESAME-style table: ``303``
+ ANEOS (in SESAME-style tables): ``4``
    + Forsterite (Stewart et al. 2019): ``400``
    + Iron (Stewart, zenodo.org/record/3866507): ``401``
    + Fe85Si15 (Stewart, zenodo.org/record/3866550): ``402``
    
The data files for the tabulated EoS can be downloaded using 
the ``examples/EoSTables/get_eos_tables.sh`` script.

To enable one or multiple EoS, the corresponding ``planetary_use_*:``
flag(s) must be set to ``1`` in the parameter file for a simulation,
along with the path to any table files, which are provided with the 
``planetary_*_table_file:`` parameters,
as detailed in :ref:`Parameters_eos` and ``examples/parameter_example.yml``.
This currently means that all EoS within each base type are prepared at once, 
which we intend to simplify in the future.

Unlike the EoS for an ideal or isothermal gas, these more complicated materials 
do not always include transformations between the internal energy, 
temperature, and entropy. At the moment, we have implemented 
\\(P(\\rho, u)\\) and \\(c_s(\\rho, u)\\), 
which is sufficient for the :ref:`planetary_sph` hydro scheme, 
but makes these materials currently incompatible with entropy-based schemes.

The Tillotson sound speed was derived using 
\\(c_s^2 = \\left. ( \\partial P / \\partial \\rho ) \\right|_S \\)
as described in 
`Kegerreis et al. (2019)  <https://doi.org/10.1093/mnras/stz1606>`_. 
Note that there is a typo in the sign of
\\(du = T dS - P dV = T dS + (P / \\rho^2) d\\rho \\) in the appendix;
the correct version was used in the derivation.
