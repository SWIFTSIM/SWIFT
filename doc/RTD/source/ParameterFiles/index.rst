.. Parameter Files
   Matthieu Schaller, 21st October 2018

Parameter Files
===============

File format and basic information
---------------------------------

The parameter file uses a format similar to the `YAML format
<https://en.wikipedia.org/wiki/YAML>`_ but reduced to only the
elements required for the SWIFT parameters. Options are given by a
name followed by a column and the value of the parameter:

.. code:: YAML

   ICs:        santa_barbara.hdf5	  
   dt_max:     1.5
   shift:      [2., 4., 5.]

Comments can be inserted anywhere and start with a hash:

.. code:: YAML

   # Descrption of the physics
   viscosity_alpha:     2.0
   dt_max:              1.5     # seconds

A typical SWIFT parameter file is split into multiple sections that
may or may not be present depending on the different configuration
options. The sections start with a label and can contain any number of
parameters:

.. code:: YAML

   Cosmology:    # Planck13
     Omega_m:        0.307
     Omega_lambda:   0.693
     Omega_b:        0.0455
     h:              0.6777
     a_begin:        0.0078125     # z = 127

The options can be integer values, floating point numbers, characters
or strings. If SWIFT expects a number and string is given, an error
will be raised. The code can also read an array of values:

.. code:: YAML

   shift:  [2., 4., 5.]
	  
Some options in the parameter file are optional and
when not provided, SWIFT will run with the default value. However, if
a compulsory parameter is missing an error will be raised at
start-up.

Finally, SWIFT outputs two YAML files at the start of a run. The first
one ``used_parameters.yml`` contains all the parameters that were used
for this run, **including all the optional parameters with their
default values**. This file can be used to start an exact copy of the
run. The second file, ``unused_parameters.yml`` contains all the
values that were not read from the parameter file. This can be used to
simplify the parameter file or check that nothing important was
ignored (for instance because the code is not configured to use some
options).

The rest of this page describes all the SWIFT parameters, split by
section. A list of all the possible parameters is kept in the file
``examples/parameter_examples.yml``.

Internal Unit System
--------------------

This section describes the units used internally by the code. This is
the system of units in which all the equations are solved. All
physical constants are converted to this system and if the ICs use a
different system (see :ref:`ICs_units_label`) the particle quantities
will be converted when read in.

The system of units is described using the value of the 5 basic units
of any system with respect to the CGS system. Instead of using a unit
of time we use a unit of velocity as this is more intuitive. Users
hence need to provide:

* a unit of length: ``UnitLength_in_cgs``,
* a unit of mass: ``UnitMass_in_cgs``,
* a unit of velocity ``UnitVelocity_in_cgs``,
* a unit of electric current ``UnitCurrent_in_cgs``,
* a unit of temperature ``UnitTemp_in_cgs``.

All these need to be expressed with respect to their cgs counter-part
(i.e. :math:`cm`, :math:`g`, :math:`cm/s`, :math:`A` and :math:`K`). Recall
that there are no h-factors in any of SWIFT's quantities; we, for instance,
use :math:`cm` and not :math:`cm/h`.

For instance to use the commonly adopted system of 10^10 Msun as a
unit for mass, mega-parsec as a unit of length and km/s as a unit of
speed, we would use:

.. code:: YAML

   # Common unit system for cosmo sims
   InternalUnitSystem:
     UnitMass_in_cgs:     1.98848e43    # 10^10 M_sun in grams
     UnitLength_in_cgs:   3.08567758e24 # 1 Mpc in centimeters
     UnitVelocity_in_cgs: 1e5           # 1 km/s in centimeters per second
     UnitCurrent_in_cgs:  1             # 1 Ampere
     UnitTemp_in_cgs:     1             # 1 Kelvin   
	  
Note that there are currently no variables in any of the SWIFT physics
schemes that make use of the unit of electric current. There is also
no incentive to use anything else than Kelvin but that makes the whole
system consistent with any possible unit system.

If one is interested in using the more humourous `FFF unit
system <https://en.wikipedia.org/wiki/FFF_system>`_ one would use

.. code:: YAML

   # FFF unit system
   InternalUnitSystem:
     UnitMass_in_cgs:     40823.3133  # 1 Firkin (fir) in grams
     UnitLength_in_cgs:   20116.8     # 1 Furlong (fur) in cm
     UnitVelocity_in_cgs: 0.01663095  # 1 Furlong (fur) per Fortnight (ftn) in cm/s
     UnitCurrent_in_cgs:  1           # 1 Ampere
     UnitTemp_in_cgs:     1           # 1 Kelvin   

The value of the physical constants in this system is left as an
exercise for the reader [#f1]_.

Cosmology
---------

When running a cosmological simulation, this section set the values of the
cosmological model. The epanded :math:`\Lambda\rm{CDM}` parameters governing the
background evolution of the Univese need to be specified here. These are:

* The reduced Hubble constant: :math:`h`: ``h``,
* The matter density parameter :math:`\Omega_m`: ``Omega_m``,
* The cosmological constant density parameter :math:`\Omega_\Lambda`: ``Omega_lambda``,
* The baryon density parameter :math:`\Omega_b`: ``Omega_b``,
* The radiation density parameter :math:`\Omega_r`: ``Omega_r``.

The last parameter can be omitted and will default to :math:`\Omega_r = 0`.

This section als specifies the start and end of the simulation expressed in
terms of scale-factors. The two parameters are:

* Initial scale-factor: ``a_begin``,
* Final scale-factor: ``a_end``.

Two additional optional parameters can be used to change the equation of
state of dark energy :math:`w(a)`. We use the evolution law :math:`w(a) =
w_0 + w_a (1 - a)`. The two parameters in the YAML file are:

* The :math:`z=0` dark energy equation of state parameter :math:`w_0`: ``w_0``
* The dark energy equation of state evolutio parameter :math:`w_a`: ``w_a``

If unspecified these parameters default to the default
:math:`\Lambda\rm{CDM}` values of :math:`w_0 = -1` and :math:`w_a = 0`.

For a Planck+13 cosmological model (ignoring radiation density as is
commonly done) and running from :math:`z=127` to :math:`z=0`, one would hence
use the following parameters:

.. code:: YAML

   Cosmology:
     a_begin:        0.0078125     # z = 127
     a_end:          1.0           # z = 0
     h:              0.6777        
     Omega_m:        0.307         
     Omega_lambda:   0.693         
     Omega_b:        0.0455        
     Omega_r:        0.            # (Optional)
     w_0:            -1.0          # (Optional)
     w_a:            0.            # (Optional)

When running a non-cosmological simulation (i.e. without the ``-c`` runtime
flag) this section of the YAML file is entirely ignored.
     
Gravity
-------

SPH
---

Time Integration
----------------

Physical Constants
------------------

For some idealised test it can be useful to overwrite the value of
some physical constants. In particular the value of the gravitational
constant. SWIFT offers an optional parameter to overwrite the value of
that constant.

.. code:: YAML

   PhysicalConstants:
     G:   1

Note that this set :math:`G` to the specified value in the internal system
of units. Setting a value of `1` when using the system of units (10^10 Msun,
Mpc, km/s) will mean that :math:`G_N=1` in these units [#f2]_ instead of the
normal value :math:`G_N=43.00927`.

This option is only used for specific tests and debugging. This entire
section of the YAML file can typically be left out.

Snapshots
---------

Some additional specific options for the snapshot outputs are described in the
following pages:

.. toctree::
   :maxdepth: 1

   output_selection

Statistics
----------

Restarts
--------

Scheduler
---------

Domain Decomposition
--------------------

.. [#f1] The thorough reader (or overly keen SWIFT tester) would find  that the speed of light is :math:`c=1.8026\times10^{12}\,\rm{fur}\,\rm{ftn}^{-1}`, Newton's contant becomes :math:`G_N=4.896735\times10^{-4}~\rm{fur}^3\,\rm{fir}^{-1}\,\rm{ftn}^{-2}` and Planck's constant turns into :math:`h=4.851453\times 10^{-34}~\rm{fur}^2\,\rm{fir}\,\rm{ftn}^{-1}`.


.. [#f2] which would translate into a constant :math:`G_N=1.5517771\times10^{-9}~cm^{3}\,g^{-1}\,s^{-2}` if expressed in the CGS system.
