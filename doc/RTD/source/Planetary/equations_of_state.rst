.. Planetary EoS
    Jacob Kegerreis, 13th March 2020

.. _planetary_eos:

Planetary Equations of State
============================
   
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
