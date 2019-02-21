.. Parameter Description
   Matthieu Schaller, 21st October 2018

.. _Parameters_basics:

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

   # Description of the physics
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

Finally, SWIFT outputs two YAML files at the start of a run. The first one
``used_parameters.yml`` contains all the parameters that were used for this run,
**including all the optional parameters left unspecified with their default
values**. This file can be used to start an exact copy of the run. The second
file, ``unused_parameters.yml`` contains all the values that were not read from
the parameter file. This can be used to simplify the parameter file or check
that nothing important was ignored (for instance because the code is not
configured to use some options).

The rest of this page describes all the SWIFT parameters, split by
section. A list of all the possible parameters is kept in the file
``examples/parameter_examples.yml``.

.. _Parameters_units:

Internal Unit System
--------------------

The ``InternalUnitSystem`` section describes the units used internally by the
code. This is the system of units in which all the equations are solved. All
physical constants are converted to this system and if the ICs use a different
system (see the snapshots' ref:`ICs_units_label` section of the documentation)
the particle quantities will be converted when read in.

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

If one is interested in using the more humorous `FFF unit
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

.. _Parameters_cosmology:

Cosmology
---------

When running a cosmological simulation, the section ``Cosmology`` sets the values of the
cosmological model. The expanded :math:`\Lambda\rm{CDM}` parameters governing the
background evolution of the Universe need to be specified here. These are:

* The reduced Hubble constant: :math:`h`: ``h``,
* The matter density parameter :math:`\Omega_m`: ``Omega_m``,
* The cosmological constant density parameter :math:`\Omega_\Lambda`: ``Omega_lambda``,
* The baryon density parameter :math:`\Omega_b`: ``Omega_b``,
* The radiation density parameter :math:`\Omega_r`: ``Omega_r``.

The last parameter can be omitted and will default to :math:`\Omega_r = 0`. Note
that SWIFT will verify on start-up that the matter content of the initial conditions
matches the cosmology specified in this section.

This section also specifies the start and end of the simulation expressed in
terms of scale-factors. The two parameters are:

* Initial scale-factor: ``a_begin``,
* Final scale-factor: ``a_end``.

Two additional optional parameters can be used to change the equation of
state of dark energy :math:`w(a)`. We use the evolution law :math:`w(a) =
w_0 + w_a (1 - a)`. The two parameters in the YAML file are:

* The :math:`z=0` dark energy equation of state parameter :math:`w_0`: ``w_0``
* The dark energy equation of state evolution parameter :math:`w_a`: ``w_a``

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

When running a non-cosmological simulation (i.e. without the ``-c`` run-time
flag) this section of the YAML file is entirely ignored.

.. _Parameters_gravity:

Gravity
-------

The behaviour of the self-gravity solver can be modified by the parameters
provided in the ``Gravity`` section. The theory document puts these parameters into the
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
:math:`\overrightarrow{a}` is the particle's acceleration. `Power et al. (2003) <http://adsabs.harvard.edu/abs/2003MNRAS.338...14P>`_ recommend using :math:`\eta=0.025`.
The last tree-related parameter is

* The tree rebuild frequency: ``rebuild_frequency``.

The tree rebuild frequency is an optional parameter defaulting to
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
* The scale below which the short-range forces are assumed to be exactly Newtonian (in units of
  the mesh cell-size multiplied by :math:`a_{\rm smooth}`) :math:`r_{\rm
  cut,min}`: ``r_cut_min`` (default: ``0.1``),

For most runs, the default values can be used. Only the number of cells along
each axis needs to be specified. The remaining three values are best described
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


.. _Parameters_SPH:

SPH
---

.. _Parameters_time_integration:

Time Integration
----------------

The ``TimeIntegration`` section is used to set some general parameters related to time
integration. In all cases, users have to provide a minimal and maximal time-step
size:

* Maximal time-step size: ``dt_max``
* Minimal time-step size: ``dt_min``

These quantities are expressed in internal units. All particles will have their
time-step limited by the maximal value on top of all the other criteria that may
apply to them (gravity acceleration, Courant condition, etc.). If a particle
demands a time-step size smaller than the minimum, SWIFT will abort with an
error message. This is a safe-guard against simulations that would never
complete due to the number of steps to run being too large.

When running a non-cosmological simulation, the user also has to provide the
time of the start and the time of the end of the simulation:

* Start time: ``time_begin``
* End time: ``time_end``

Both are expressed in internal units. The start time is typically set to ``0``
but SWIFT can handle any value here. For cosmological runs, these values are
ignored and the start- and end-points of the runs are specified by the start and
end scale-factors in the cosmology section of the parameter file.

Additionally, when running a cosmological volume, advanced users can specify the
value of the dimensionless pre-factor entering the time-step condition linked
with the motion of particles with respect to the background expansion and mesh
size. See the theory document for the exact equations.

* Dimensionless pre-factor of the maximal allowed displacement:
  ``max_dt_RMS_factor`` (default: ``0.25``)

This value rarely needs altering.

A full time-step section for a non-cosmological run would be:

.. code:: YAML

  TimeIntegration:
    time_begin:   0    # Start time in internal units.
    time_end:     10.  # End time in internal units.
    dt_max:       1e-2
    dt_min:       1e-6

Whilst for a cosmological run, one would need:

.. code:: YAML

  TimeIntegration:
    dt_max:            1e-4
    dt_min:            1e-10
    max_dt_RMS_factor: 0.25     # Default optional value

.. _Parameters_ICs:

Initial Conditions
------------------

The ``InitialConditions`` section of the parameter file contains all the options related to
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


.. _Parameters_constants:

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

.. _Parameters_snapshots:

Snapshots
---------

The ``Snapshots`` section of the parameter file contains all the options related to
the dump of simulation outputs in the form of HDF5 :ref:`snapshots`. The main
parameter is the base name that will be used for all the outputs in the run:

* The base name of the HDF5 snapshots: ``basename``.

This name will then be appended by an under-score and 4 digits followed by
``.hdf5`` (e.g. ``base_name_1234.hdf5``). The 4 digits are used to label the
different outputs, starting at ``0000``. In the default setup the digits simply
increase by one for each snapshot. However, if the optional parameter
``int_time_label_on`` is switched on, then we use 6 digits and these will the
physical time of the simulation rounded to the nearest integer
(e.g. ``base_name_001234.hdf5``) [#f3]_.

The time of the first snapshot is controlled by the two following options:

* Time of the first snapshot (non-cosmological runs): ``time_first``,
* Scale-factor of the first snapshot (cosmological runs): ``scale_factor_first``.

One of those two parameters has to be provided depending on the type of run. In
the case of non-cosmological runs, the time of the first snapshot is expressed
in the internal units of time. Users also have to provide the difference in time
(or scale-factor) between consecutive outputs:

* Time difference between consecutive outputs: ``delta_time``.

In non-cosmological runs this is also expressed in internal units. For
cosmological runs, this value is *multiplied* to obtain the
scale-factor of the next snapshot. This implies that the outputs are
equally space in :math:`\log(a)` (See :ref:`Output_list_label` to have
snapshots not regularly spaced in time).

When running the code with structure finding activated, it is often
useful to have a structure catalog written at the same simulation time
as the snapshots. To activate this, the following parameter can be
switched on:

* Run VELOCIraptor every time a snapshot is dumped: ``invoke_stf``
  (default: ``0``).

This produces catalogs using the options specified for the stand-alone
VELOCIraptor outputs (see the section :ref:`Parameters_structure_finding`) but
with a base name and output number that matches the snapshot name
(e.g. ``stf_base_name_1234.hdf5``) irrespective of the name specified in the
section dedicated to VELOCIraptor. Note that the invocation of VELOCIraptor at
every dump is done additionally to the stand-alone dumps that can be specified
in the corresponding section of the YAML parameter file.

Users can optionally specify the level of compression used by the HDF5 library
using the parameter:

* GZIP compression level of the HDF5 arrays: ``compression`` (default: ``0``).

The default level of ``0`` implies no compression and values have to be in the
range :math:`[0-9]`. This integer is passed to the i/o library and used for the
loss-less GZIP compression algorithm. Higher values imply higher compression but
also more time spent deflating and inflating the data. Note that up until HDF5
1.10.x this option is not available when using the MPI-parallel version of the
i/o routines.

Finally, it is possible to specify a different system of units for the snapshots
than the one that was used internally by SWIFT. The format is identical to the
one described above (See the :ref:`Parameters_units` section) and read:

* a unit of length: ``UnitLength_in_cgs`` (default: ``InternalUnitSystem:UnitLength_in_cgs``),
* a unit of mass: ``UnitMass_in_cgs`` (default: ``InternalUnitSystem:UnitMass_in_cgs``),
* a unit of velocity ``UnitVelocity_in_cgs`` (default: ``InternalUnitSystem:UnitVelocity_in_cgs``),
* a unit of electric current ``UnitCurrent_in_cgs`` (default: ``InternalUnitSystem:UnitCurrent_in_cgs``),
* a unit of temperature ``UnitTemp_in_cgs`` (default: ``InternalUnitSystem:UnitTemp_in_cgs``).

When un-specified, these all take the same value as assumed by the internal
system of units. These are rarely used but can offer a practical alternative to
converting data in the post-processing of the simulations.

For a standard cosmological run with structure finding activated, the
full section would be:

.. code:: YAML

   Snapshots:
     basename:            output
     scale_factor_first:  0.02    # z = 49
     delta_time:          1.02
     invoke_stf:          1

Showing all the parameters for a basic hydro test-case, one would have:

.. code:: YAML

   Snapshots:
     basename:            sedov
     time_first:          0.01
     delta_time:          0.005
     invoke_stf:          0
     int_time_label_on:   0
     compression:         3
     UnitLength_in_cgs:   1.  # Use cm in outputs
     UnitMass_in_cgs:     1.  # Use grams in outputs
     UnitVelocity_in_cgs: 1.  # Use cm/s in outputs
     UnitCurrent_in_cgs:  1.  # Use Ampere in outputs
     UnitTemp_in_cgs:     1.  # Use Kelvin in outputs

Some additional specific options for the snapshot outputs are described in the
following pages:

* :ref:`Output_list_label` (to have snapshots not evenly spaced in time),
* :ref:`Output_selection_label` (to select what particle fields to write).


.. _Parameters_statistics:

Statistics
----------

Some additional specific options for the statistics outputs are described in the
following page:

* :ref:`Output_list_label` (to have statistics outputs not evenly spaced in time).

.. _Parameters_restarts:

Restarts
--------

SWIFT can write check-pointing files and restart from them. The behaviour of
this mechanism is driven by the options in the ``Restarts`` section of the YAML
parameter file. All the parameters are optional but default to values that
ensure a reasonable behaviour.

* Whether or not to enable the dump of restart files: ``enable`` (default:
  ``1``).

This parameter acts a master-switch for the check-pointing capabilities. All the
other options require the ``enable`` parameter to be set to ``1``.

* Whether or not to save a copy of the previous set of check-pointing files:
  ``save`` (default: ``1``),
* Whether or not to dump a set of restart file on regular exit: ``onexit``
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
here. The files themselves are named ``basename_000001.rst`` where the basename
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

.. _Parameters_scheduler:

Scheduler
---------

The Scheduler section contains various parameters that control how the cell
tree is configured and defines some values for the related tasks.  In general
these should be considered as tuning parameters, both for speed and memory
use.

.. code:: YAML

   nr_queues: 0

Defines the number of task queues used. These are normally set to one per
thread and should be at least that number.

A number of parameters decide how the cell tree will be split into sub-cells,
according to the number of particles and their expected interaction count,
and the type of interaction. These are:

.. code:: YAML

  cell_max_size:             8000000
  cell_sub_size_pair_hydro:  256000000
  cell_sub_size_self_hydro:  32000
  cell_sub_size_pair_grav:   256000000
  cell_sub_size_self_grav:   32000
  cell_sub_size_pair_stars:  256000000
  cell_sub_size_self_stars:  32000
  cell_split_size:           400

when possible cells that exceed these constraints will be split into a further
level of sub-cells. So for instance a sub-cell should not contain more than
400 particles (this number defines the scale of most `N*N` interactions).

To control the number of self-gravity tasks we have the parameter:

.. code:: YAML

  cell_subdepth_diff_grav:   4

which stops these from being done at the scale of the leaf cells, of which
there can be a large number. In this case cells with gravity tasks must be at
least 4 levels above the leaf cells (when possible).

Extra space is required when particles are created in the system (to the time
of the next rebuild). These are controlled by:

.. code:: YAML

  cell_extra_parts:          0
  cell_extra_gparts:         0
  cell_extra_sparts:         400


The number of top-level cells is controlled by the parameter:

.. code:: YAML

  max_top_level_cells:       12

this is the number per dimension, we will have 12x12x12 cells. There must be
at least 3 top-level cells per dimension.

The number of top-level cells should be set so that the number of particles
per cell is not too large, this is particularly important when using MPI as
this defines the maximum size of cell exchange and also the size of non-local
cells (these are used for cell interactions with local cells), which can have
a large influence on memory use. Best advice for this is to at least scale for
additional nodes.

The memory used for holding the task and task-link lists needs to be
pre-allocated, but cannot be pre-calculated, so we have the two parameters:

.. code:: YAML

  tasks_per_cell:            0.0
  links_per_tasks:           10

which are guesses at the mean numbers of tasks per cell and number of links
per task. The tasks_per_cell value will be conservatively guessed when set to
0.0, but you will be able to save memory by setting a value. The way to get a
better estimate is to run SWIFT with verbose reporting on (```--verbose=1```)
and check for the lines that report the ```per cell``` or with MPI ``maximum
per cell``` values. This number can vary as the balance between MPI ranks
does, so it is probably best to leave some head room.

If these are exceeded you should get an obvious error message.

Finally the parameter:

.. code:: YAML

  mpi_message_limit:         4096

Defines the size (in bytes) below which MPI communication will be sent using
non-buffered calls. These should have lower latency, but how that works or
is honoured is an implementation question.


.. _Parameters_domain_decomposition:

Domain Decomposition:
---------------------

This section determines how the top-level cells are distributed between the
ranks of an MPI run. An ideal decomposition should result in each rank having
a similar amount of work to do, so that all the ranks complete at the same
time. Achieving a good balance requires that SWIFT is compiled with either the
ParMETIS or METIS libraries. ParMETIS is an MPI version of METIS, so is
preferred for performance reasons.

When we use ParMETIS/METIS the top-level cells of the volume are considered as
a graph, with a cell at each vertex and edges that connect the vertices to all
the neighbouring cells (so we have 26 edges connected to each vertex).
Decomposing such a graph into domains is known as partitioning, so in SWIFT we
refer to domain decomposition as partitioning.

This graph of cells can have weights associated with the vertices and the
edges. These weights are then used to guide the partitioning, seeking to
balance the total weight of the vertices and minimize the weights of the edges
that are cut by the domain boundaries (known as the edgecut). We can consider
the edge weights as a proxy for the exchange of data between cells, so
minimizing this reduces communication.

The Initial Partition:
^^^^^^^^^^^^^^^^^^^^^^

When SWIFT first starts it reads the initial conditions and then does an
initial distribution of the top-level cells. At this time the only information
available is the cell structure and, by geometry, the particles each cell
should contain. The type of partitioning attempted is controlled by the::

  DomainDecomposition:
    initial_type:

parameter. Which can have the values *memory*, *region*, *grid* or
*vectorized*:


    * *memory*

    This is the default if METIS or ParMETIS is available. It performs a
    partition based on the memory use of all the particles in each cell,
    attempting to equalize the memory used by all the ranks.
    How successful this attempt is depends on the granularity of cells and particles
    and the number of ranks, clearly if most of the particles are in one cell,
    or a small region of the volume, balance is impossible or
    difficult. Having more top-level cells makes it easier to calculate a
    good distribution (but this comes at the cost of greater overheads).

    * *region*

    The one other METIS/ParMETIS option is "region". This attempts to assign equal
    numbers of cells to each rank, with the surface area of the regions minimised
    (so we get blobs, rather than rectangular volumes of cells).

If ParMETIS and METIS are not available two other options are possible, but
will give a poorer partition:

    * *grid*

    Split the cells into a number of axis aligned regions. The number of
    splits per axis is controlled by the::

       initial_grid

    parameter. It takes an array of three values. The product of these values
    must equal the number of MPI ranks. If not set a suitable default will be used.

    * *vectorized*

    Allocate the cells on the basis of proximity to a set of seed
    positions. The seed positions are picked every nranks along a vectorized
    cell list (1D representation). This is guaranteed to give an initial
    partition for all cases when the number of cells is greater equal to the
    number of MPI ranks, so can be used if the others fail. Don't use this.

If ParMETIS and METIS are not available then only an initial partition will be
performed. So the balance will be compromised by the quality of the initial
partition.

Repartitioning:
^^^^^^^^^^^^^^^

When ParMETIS or METIS is available we can consider adjusting the balance
during the run, so we can improve from the initial partition and also track
changes in the run that require a different balance. The initial partition is
usually not optimal as although it may have balanced the distribution of
particles it has not taken account of the fact that different particles types
require differing amounts of processing and we have not considered that we
also need to do work requiring communication between cells. This latter point
is important as we are running an MPI job, as inter-cell communication may be
very expensive.

There are a number of possible repartition strategies which are defined using
the::

  DomainDecomposition:
    repartition_type:

parameter. The possible values for this are *none*, *fullcosts*, *edgecosts*,
*memory*, *timecosts*.

    * *none*

    Rather obviously, don't repartition. You are happy to run with the
    initial partition.

    * *fullcosts*

    Use computation weights derived from the running tasks for the vertex and
    edge weights. This is the default.

    * *edgecosts*

    Only use computation weights derived from the running tasks for the edge
    weights.

    * *memory*

    Repeat the initial partition with the current particle positions
    re-balancing the memory use.

    * *timecosts*

    Only use computation weights derived from the running tasks for the vertex
    weights and the expected time the particles will interact in the cells as
    the edge weights. Using time as the edge weight has the effect of keeping
    very active cells on single MPI ranks, so can reduce MPI communication.

The computation weights are actually the measured times, in CPU ticks, that
tasks associated with a cell take. So these automatically reflect the relative
cost of the different task types (SPH, self-gravity etc.), and other factors
like how well they run on the current hardware and are optimized by the
compiler used, but this means that we have a constraint on how often we can
consider repartitioning, namely when all (or nearly all) the tasks of the
system have been invoked in a step. To control this we have the::

    minfrac:     0.9

parameter. Which defines the minimum fraction of all the particles in the
simulation that must have been actively updated in the last step, before
repartitioning is considered.

That then leaves the question of when a run is considered to be out of balance
and should benefit from a repartition. That is controlled by the::

    trigger:          0.05

parameter. This value is the CPU time difference between MPI ranks, as a
fraction, if less than this value a repartition will not be
done. Repartitioning can be expensive not just in CPU time, but also because
large numbers of particles can be exchanged between MPI ranks, so is best
avoided.

If you are using ParMETIS there additional ways that you can tune the
repartition process.

METIS only offers the ability to create a partition from a graph, which means
that each solution is independent of those that have already been made, that
can make the exchange of particles very large (although SWIFT attempts to
minimize this), however, using ParMETIS we can use the existing partition to
inform the new partition, this has two algorithms that are controlled using::

    adaptive:         1

which means use adaptive repartition, otherwise simple refinement. The
adaptive algorithm is further controlled by the::

    itr:              100

parameter, which defines the ratio of inter node communication time to data
redistribution time, in the range 0.00001 to 10000000.0. Lower values give
less data movement during redistributions. The best choice for these can only
be determined by experimentation (the gains are usually small, so not really
recommended).

Finally we have the parameter::

    usemetis:         0

Forces the use of the METIS API, probably only useful for developers.

**Fixed cost repartitioning:**

So far we have assumed that repartitioning will only happen after a step that
meets the `minfrac:` and `trigger:` criteria, but we may want to repartition
at some arbitrary steps, and indeed do better than the initial partition
earlier in the run. This can be done using *fixed cost* repartitioning.

Fixed costs are output during each repartition step into the file
`partition_fixed_costs.h`, this should be created by a test run of your your
full simulation (with possibly with a smaller volume, but all the physics
enabled). This file can then be used to replace the same file found in the
`src/` directory and SWIFT should then be recompiled. Once you have that, you
can use the parameter::

    use_fixed_costs:  1

to control whether they are used or not. If enabled these will be used to
repartition after the second step, which will generally give as good a
repartition immediately as you get at the first unforced repartition.

Also once these have been enabled you can change the `trigger:` value to
numbers greater than 2, and repartitioning will be forced every `trigger`
steps. This latter option is probably only useful for developers, but tuning
the second step to use fixed costs can give some improvements.

.. _Parameters_structure_finding:

Structure finding (VELOCIraptor)
--------------------------------


.. [#f1] The thorough reader (or overly keen SWIFT tester) would find  that the speed of light is :math:`c=1.8026\times10^{12}\,\rm{fur}\,\rm{ftn}^{-1}`, Newton's constant becomes :math:`G_N=4.896735\times10^{-4}~\rm{fur}^3\,\rm{fir}^{-1}\,\rm{ftn}^{-2}` and Planck's constant turns into :math:`h=4.851453\times 10^{-34}~\rm{fur}^2\,\rm{fir}\,\rm{ftn}^{-1}`.


.. [#f2] which would translate into a constant :math:`G_N=1.5517771\times10^{-9}~cm^{3}\,g^{-1}\,s^{-2}` if expressed in the CGS system.

.. [#f3] This feature only makes sense for non-cosmological runs for which the
         internal time unit is such that when rounded to the nearest integer a
	 sensible number is obtained. A use-case for this feature would be to
	 compare runs over the same physical time but with different numbers of
	 snapshots. Snapshots at a given time would always have the same set of
	 digits irrespective of the number of snapshots produced before.

