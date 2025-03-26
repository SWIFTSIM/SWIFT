.. Planetary EoS
    Jacob Kegerreis, 14th July 2022

.. _planetary_eos:

Planetary Equations of State
============================

Configuring SWIFT with the ``--with-equation-of-state=planetary`` and
``--with-hydro=planetary`` options enables the use of multiple
equations of state (EoS).
Every SPH particle then requires and carries the additional ``MaterialID`` flag
from the initial conditions file. This flag indicates the particle's material
and which EoS it should use.

If you have another EoS that you would like us to add, then just let us know!

It is important to check that the EoS you use are appropriate
for the conditions in the simulation that you run.
Please follow the original sources of these EoS for more information and
to check the regions of validity. If an EoS sets particles to have a pressure
of zero, then particles may end up overlapping, especially if the gravitational
softening is very small.

So far, we have implemented several Tillotson, ANEOS, SESAME,
and Hubbard \& MacFarlane (1980) materials, with more on the way.
Custom materials in SESAME-style tables can also be provided.
The material's ID is set by a somewhat arbitrary base type ID
(multiplied by 100) plus an individual value, matching our code for making
planetary initial conditions, `WoMa  <https://github.com/srbonilla/WoMa>`_:

+ Ideal gas: ``0``
    + Default (Set :math:`\gamma` using ``--with-adiabatic-index``, default 5/3): ``0``
+ Tillotson (Melosh, 2007): ``1``
    + Iron: ``100``
    + Granite: ``101``
    + Water: ``102``
    + Basalt: ``103``
+ Hubbard \& MacFarlane (1980): ``2``
    + Hydrogen-helium atmosphere: ``200``
    + Ice H20-CH4-NH3 mix: ``201``
    + Rock SiO2-MgO-FeS-FeO mix: ``202``
+ SESAME (and others in similar-style tables): ``3``
    + Iron (2140): ``300``
    + Basalt (7530): ``301``
    + Water (7154): ``302``
    + Senft \& Stewart (2008) water: ``303``
+ ANEOS (in SESAME-style tables): ``4``
    + Forsterite (Stewart et al. 2019): ``400``
    + Iron (Stewart, zenodo.org/record/3866507): ``401``
    + Fe85Si15 (Stewart, zenodo.org/record/3866550): ``402``
+ Custom (in SESAME-style tables): ``9``
    + User-provided custom material(s): ``900``, ``901``, ..., ``909``

The data files for the tabulated EoS can be downloaded using
the ``examples/Planetary/EoSTables/get_eos_tables.sh`` script.

To enable one or multiple EoS, the corresponding ``planetary_use_*:``
flag(s) must be set to ``1`` in the parameter file for a simulation,
along with the path to any table files, which are set by the
``planetary_*_table_file:`` parameters,
as detailed in :ref:`Parameters_eos` and ``examples/parameter_example.yml``.

Unlike the EoS for an ideal or isothermal gas, these more complicated materials
do not always include transformations between the internal energy,
temperature, and entropy. At the moment, we have implemented
:math:`P(\rho, u)` and :math:`c_s(\rho, u)` (and more in some cases),
which is sufficient for the :ref:`planetary_sph` hydro scheme,
but some materials may thus currently be incompatible with
e.g. entropy-based schemes.

The Tillotson sound speed was derived using
:math:`c_s^2 = \left. ( \partial P / \partial \rho ) \right|_S`
as described in
`Kegerreis et al. (2019)  <https://doi.org/10.1093/mnras/stz1606>`_.
Note that there is a typo in the sign of
:math:`du = T dS - P dV = T dS + (P / \rho^2) d\rho` in the appendix,
but the correct version was used in the actual derivation.

The ideal gas uses the same equations detailed in :ref:`equation_of_state`.

The data files for the tabulated EoS can be downloaded using
the ``examples/EoSTables/get_eos_tables.sh`` script.

The format of the data files for SESAME, ANEOS, and similar-EoS tables
is similar to the SESAME 301 (etc) style. The file contents are:

.. code-block:: python

    # header (12 lines)
    version_date                                                (YYYYMMDD)
    num_rho  num_T
    rho[0]   rho[1]  ...  rho[num_rho]                          (kg/m^3)
    T[0]     T[1]    ...  T[num_T]                              (K)
    u[0, 0]                 P[0, 0]     c[0, 0]     s[0, 0]     (J/kg, Pa, m/s, J/K/kg)
    u[1, 0]                 ...         ...         ...
    ...                     ...         ...         ...
    u[num_rho-1, 0]         ...         ...         ...
    u[0, 1]                 ...         ...         ...
    ...                     ...         ...         ...
    u[num_rho-1, num_T-1]   ...         ...         s[num_rho-1, num_T-1]

The ``version_date`` must match the value in the ``sesame.h`` ``SESAME_params``
objects, so we can ensure that any version updates work with the git repository.
This is ignored for custom materials.
The header contains a first line that gives the material name, followed by the
same 11 lines printed here to describe the contents.
