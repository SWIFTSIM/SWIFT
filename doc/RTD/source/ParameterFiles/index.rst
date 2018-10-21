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
section.

InternalUnitSystem
------------------

PhysicalConstants
-----------------

Cosmology
---------

Gravity
-------

SPH
---

TimeIntegration
---------------

Snapshots
---------

Statistics
----------

Restarts
--------

Scheduler
---------

DomainDecomposition
-------------------


