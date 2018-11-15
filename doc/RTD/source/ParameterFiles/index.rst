.. Parameter Files
   Matthieu Schaller, 21st October 2018

.. _Parameter_File_label:

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

The last parameter can be omitted and will default to :math:`\Omega_r = 0`. Note
that SWIFT will verify on start-up that the matter content of the initial conditions
matches the cosmology specified in this section.

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

The behaviour of the self-gravity solver can be modifed by the parameters
provided in this section. The theory document puts these parameters into the
context of the equations being solved. We give a brief overview here.

* The Plummer-equivalent co-moving softening length used for all particles :math:`\epsilon_{com}`: ``comoving_softening``,
* The Plummer-equivalent maximal physical softening length used for all particles :math:`\epsilon_{max}`: ``comoving_softening``, 

At any redshift :math:`z`, the Plummer-equivalent softening length used by the
code will be :math:`\epsilon=\min(\epsilon_{max},
\frac{\epsilon_{com}}{z+1})`. This is expressed in internal units.

* The opening angle (multipole acceptance criterion) used in the FMM :math:`\theta`: ``theta``,
* The time-step size pre-factor :math:`\eta`: ``eta``,
  
The time-step of a given particle is given by :math:`\Delta t =
\eta\sqrt{\frac{\epsilon}{|\overrightarrow{a}|}}`, where
:math:`\overrightarrow{a}` is the particle's acceleration. Power et al. (2003) recommend using :math:`\eta=0.025`.
The last tree-related parameter is

* The tree rebuild frequency: ``rebuild_frequency``.

Thqe tree rebuild frequency is an optional parameter defaulting to
:math:`0.01`. It is used to trigger the re-construction of the tree every time a
fraction of the particles have been integrated (kicked) forward in time.

Simulations using periodic boundary conditions use additional parameters for the
Particle-Mesh part of the calculation. The last three are optional:

* The number cells along each axis of the mesh :math:`N`: ``mesh_side_length``,
* The mesh smoothing scale in units of the mesh cell-size :math:`a_{\rm
  smooth}`: ``a_smooth`` (default: ``1.25``),
* The scale above which the short-range forces are assumed to be 0 (in units of
  the mesh cell-size multiplied by :math:`a_{\rm smooth}`) :math:`r_{\rm
  cut,max}`: ``r_cut_max`` (default: ``4.5``),
* The scale bewlo which the short-range forces are assumed to be exactly Newtonian (in units of
  the mesh cell-size multiplied by :math:`a_{\rm smooth}`) :math:`r_{\rm
  cut,min}`: ``r_cut_min`` (default: ``0.1``),
  
For most runs, the default values can be used. Only the number of cells along
each axis needs to be sepcified. The remaining three values are best described
in the context of the full set of equations in the theory documents.
  
As a summary, here are the values used for the EAGLE :math:`100^3~{\rm Mpc}^3`
simulation:

.. code:: YAML
	  
   # Parameters for the self-gravity scheme for the EAGLE-100 box
   Gravity:
     eta:          0.025              
     theta:        0.7                
     comoving_softening:     0.0026994  # 0.7 proper kpc at z=2.8.
     max_physical_softening: 0.0007     # 0.7 proper kpc
     rebuild_frequency:      0.01       # Default optional value
     mesh_side_length:       512       
     a_smooth:     1.25                 # Default optional value
     r_cut_max:    4.5                  # Default optional value
     r_cut_min:    0.1                  # Default optional value

      
SPH
---

Time Integration
----------------

Initial Conditions
------------------

This section of the parameter file contains all the options related to
the initial conditions. The main two parameters are

* The name of the initial conditions file: ``file_name``,
* Whether the problem uses periodic boundary conditions or not: ``periodic``.

The file path is relative to where the code is being executed. These
parameters can be complemented by some optional values to drive some
specific behaviour of the code.

* Whether to generate gas particles from the DM particles: ``generate_gas_in_ics`` (default: ``0``),
* Whether to activate an additional clean-up of the SPH smoothing lengths: ``cleanup_smoothing_lengths`` (default: ``0``)

The procedure used to generate gas particles from the DM ones is
outlined in the theory documents and is too long for a full
description here.  The cleaning of the smoothing lengths is an
expensive operation but can be necessary in the cases where the
initial conditions are of poor quality and the values of the smoothing
lengths are far from the values they should have.

When starting from initial conditions created for Gadget, some
additional flags can be used to convert the values from h-full to
h-free and remove the additional :math:`\sqrt{a}` in the velocities:

* Whether to re-scale all the fields to remove powers of h from the quantities: ``cleanup_h_factors`` (default: ``0``),
* Whether to re-scale the velocities to remove the :math:`\sqrt{a}` assumed by Gadget : ``cleanup_velocity_factors`` (default: ``0``).

The h-factors are self-consistently removed according to their units
and this is applied to all the quantities irrespective of particle
types. The correct power of ``h`` is always calculated for each
quantity.

Finally, SWIFT also offers these options:

* A factor to re-scale all the smoothing-lengths by a fixed amount: ``smoothing_length_scaling`` (default: ``1.``),
* A shift to apply to all the particles: ``shift`` (default: ``[0.0,0.0,0.0]``),
* Whether to replicate the box along each axis: ``replicate`` (default: ``1``).

The shift is expressed in internal units. The option to replicate the
box is especially useful for weak-scaling tests. When set to an
integer >1, the box size is multiplied by this integer along each axis
and the particles are duplicated and shifted such as to create exact
copies of the simulation volume.

The full section to start a DM+hydro run from Gadget DM-only ICs would
be:

.. code:: YAML

   InitialConditions:
     file_name:  my_ics.hdf5
     periodic:                    1
     cleanup_h_factors:           1
     cleanup_velocity_factors:    1
     generate_gas_in_ics:         1  
     cleanup_smoothing_lengths:   1  

  
Physical Constants
------------------

For some idealised test it can be useful to overwrite the value of
some physical constants; in particular the value of the gravitational
constant. SWIFT offers an optional parameter to overwrite the value of
:math:`G_N`. 

.. code:: YAML

   PhysicalConstants:
     G:   1

Note that this set :math:`G` to the specified value in the internal system
of units. Setting a value of `1` when using the system of units (10^10 Msun,
Mpc, km/s) will mean that :math:`G_N=1` in these units [#f2]_ instead of the
normal value :math:`G_N=43.00927`.

This option is only used for specific tests and debugging. This entire
section of the YAML file can typically be left out. More constants may
be handled in the same way in future versions.

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

SWIFT can write check-pointing files and restart from them. The behaviour of
this mechanism is driven by the options int he `Restarts` section of the YAML
parameter file. All the parameters are optional but default to values that
ensure a reasonable behaviour. 

* Wether or not to enable the dump of restart files: ``enable`` (default:
  ``1``).

This parameter acts a master-switch for the check-pointing capabilities. All the
other options require the ``enable`` parameter to be set to ``1``.
  
* Wether or not to save a copy of the previous set of check-pointing files:
  ``save`` (default: ``1``),
* Wether or not to dump a set of restart file on regular exit: ``onexit``
  (default: ``0``),
* The wall-clock time in hours between two sets of restart files:
  ``delta_hours`` (default: ``6.0``).
  
Note that there is no buffer time added to the ``delta_hours`` value. If the
system's batch queue run time limit is set to 6 hours, the user must specify a
smaller value to allow for enough time to safely dump the check-point files.

* The sub-directory in which to store the restart files: ``subdir`` (default:
  ``restart``),
* The basename of the restart files: ``basename`` (default: ``swift``)

If the directory does not exist, SWIFT will create it.  When resuming a run,
SWIFT, will look for files with the name provided in the sub-directory specified
here. The files themselves are named ``basename_000001.rst`` where the basenme
is replaced by the user-specified name and the 6-digits number corresponds to
the MPI-rank. SWIFT writes one file per MPI rank. If the ``save`` option has
been activated, the previous set of restart files will be named
``basename_000000.rst.prev``.

SWIFT can also be stopped by creating an empty file called ``stop`` in the
directory where the code runs. This will make SWIFT dump a fresh set of restart
file (irrespective of the specified ``delta_time`` between dumps) and exit
cleanly. One parameter governs this behaviour:

* Number of steps between two checks for the presence of a ``stop`` file:
  ``stop_steps`` (default: ``100``).

The default value is chosen such that SWIFT does not need to poll the
file-system to often, which can take a significant amount of time on distributed
systems. For runs where the small time-steps take a much larger amount of time,
a smaller value is recommended to allow for a finer control over when the code
can be stopped.

Finally, SWIFT can automatically stop after a specified amount of wall-clock
time. The code can also run a command when exiting in this fashion, which can be
used, for instance, to interact with the batch queue system:

* Maximal wall-clock run time in hours: ``max_run_time`` (default: ``24.0``),
* Whether or not to run a command on exit: ``resubmit_on_exit`` (default:
  ``0``),
* The command to run on exit: ``resubmit_command`` (default: ``./resub.sh``).

Note that no check is performed on the validity of the command to run. SWIFT
simply calls ``system()`` with the user-specified command.

To run SWIFT, dumping check-pointing files every 6 hours and running for 24
hours after which a shell command will be run, one would use:

.. code:: YAML
	  
  Restarts:
    enable:             1          
    save:               1          # Keep copies
    onexit:             0          
    subdir:             restart    # Sub-directory of the directory where SWIFT is run
    basename:           swift      
    delta_hours:        6.0        
    stop_steps:         100        
    max_run_time:       24.0       # In hours 
    resubmit_on_exit:   1          
    resubmit_command:   ./resub.sh 



Scheduler
---------

Domain Decomposition
--------------------

.. [#f1] The thorough reader (or overly keen SWIFT tester) would find  that the speed of light is :math:`c=1.8026\times10^{12}\,\rm{fur}\,\rm{ftn}^{-1}`, Newton's contant becomes :math:`G_N=4.896735\times10^{-4}~\rm{fur}^3\,\rm{fir}^{-1}\,\rm{ftn}^{-2}` and Planck's constant turns into :math:`h=4.851453\times 10^{-34}~\rm{fur}^2\,\rm{fir}\,\rm{ftn}^{-1}`.


.. [#f2] which would translate into a constant :math:`G_N=1.5517771\times10^{-9}~cm^{3}\,g^{-1}\,s^{-2}` if expressed in the CGS system.
