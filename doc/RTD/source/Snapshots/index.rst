.. Snapshots
   Matthieu Schaller, 5th January 2019

.. _snapshots:

Snapshots
=========

The snapshots are stored using the HDF5 format and are almost compatible with
Gadget-2 (fully compatible outside of cosmological runs). They do, however,
contain a large set of extensions including units, meta-data about the code and
runs as well as facilities to quickly access the particles in a specific region
of the simulation volume.

Header
------

The header group (``/Header``) contains basic information about the simulation
and the number of particles. For compatibility reasons, the SWIFT snapshot's
header format is identical to the Gadget-2 one with small additions.

In addition to the standard quantities, the header contains a field ``Code``
that is always set to the string ``SWIFT``, which can be used to identify
SWIFT-generated snapshots and hence make use of all the extensions to the file
format described below.

The most important quantity of the header is the array ``NumPart_Total`` which
contains the number of particles of each type in this snapshot. This is an array
of 6 numbers; one for each of the supported types. Following the Gadget-2
convention, if that number is larger than 2^31, SWIFT will use the
``NumPart_HighWord`` field to store the high-word bits of the total number of
particles. The field ``NumPart_ThisFile`` contains the number of particles in
this sub-snapshot file when the user asked for distributed snapshots (see
:ref:`Parameters_snapshots`); otherwise it contains the same information as
``NumPart_Total``. Note however, that there is no high-word for this field. We
store it as a 64-bits integer [#f1]_. The field ``NumFilesPerSnapshot`` specifies the
number of sub-snapshot files (always 1 unless a distributed snapshot was asked
for) and ``ThisFile`` the id of that specific file (always 0 unless a distributed
snapshot was asked for). 

The field ``TotalNumberOfParticles`` gives the total number of particles of each type
as a 64 bit integer. This allows the total number of particles to be read directly
with no calculation required even if there are 2^31 or more particles. This field is
equal to ``NumPart_ThisFile`` if the snapshot is not distributed over multiple files.

The field ``InitialMassTable`` contains the *mean* initial mass of each of the
particle types present in the initial conditions. This can be used as estimator
of the mass resolution of the run. The masses are expressed in internal units.

The field ``OutputType`` contains information about the kind of output this
snapshot is. The possible values are:

+---------------------+-----------------------------------------------------+
| OutputType          | Definition                                          |
+=====================+=====================================================+
| ``FullVolume``      | Regular vanilla snapshot                            |
+---------------------+-----------------------------------------------------+
| ``SubSampled``      | Snapshot where some particle types were sub-sampled |
+---------------------+-----------------------------------------------------+
| ``LineOfSight``     | Line-of-sight snapshot                              |
+---------------------+-----------------------------------------------------+
| ``FOF``             | Friends-Of-Friends Halo Catalogue                   |
+---------------------+-----------------------------------------------------+


The ``RunName`` field contains the name of the simulation that was specified as
the ``run_name`` in the :ref:`Parameters_meta_data` section of the YAML
parameter file.

The ``System`` field contains the name of the machine where the MPI rank 0 was
placed. This name is whatever UNIX's ``gethostname()`` function returns on that
system. Similarly, the ``SnapshotDate`` field contains the date and time when
the file was written.

The ``TimeBase_dloga`` field contains the change in logarithm of the
scale-factor corresponding to a time-step of length 1 on the integer
time-line. This is the smallest time-step size that the code can use. This field
is zero in non-cosmological runs. Similarly, the field ``TimeBase_dt`` contains
the smallest time-step size (in internal units) that the code can take. This
would be the increase in time a particle in the time-bin one would have. Note
that in cosmological runs this quantity evolves with redshift as the (logarithm
of the) scale-factor is used on the integer time-line.

The field ``SelectOutput`` will contain the name of the
:ref:`Output_selection_label` used for this specific output and will take the value
``Default`` if no such selection (or the default one) was used.

If a sub-sampling of the particle fields was used, then the header additionally
contains a field describing the fraction of the particles of each type that were
written to the snapshot. Note, however, that when sub-sampling the fields 
``NumPart_Total``, ``NumPart_HighWord``, and ``NumPart_ThisFile`` contain the number
of particles actually written (i.e. after sub-sampling), not the total number of
particles in the run.

The field ``CanHaveTypes`` contains information about whether a given particle
type is to be expected in snapshots of the run. For instance, a simulation with
star formation switched on, the code may not have formed a star yet but might in
future snapshots. This allows reading tools to distinguish fields they will
never expect to find in a given simulation from fields that may be present in
other outputs.

Finally, the shift that may have been applied to all the particles upon reading
the ICs (See :ref:`Parameters_ICs`) is added to the header in the field
``Shift``. This is expressed in internal units.

Meta-data about the code and run
--------------------------------

Several groups at the root of the files only contain attributes and are used to
store some meta-data about the simulation and the code itself.

Code
~~~~

The group ``/Code`` contains basic information about the version of the code
that was used to run the simulation that dumped this snapshot. Versions of the
libraries used to compile the code as well as information about the compiler and
the flags used are stored. The most important element here are the git SHA and
configuration parameters of the code. Alongside the compiler flags, policies and
used parameters, these allow to reproduce exactly an older run.

Cosmology
~~~~~~~~~

The group ``/Cosmology`` contains information about the cosmological model used
for this simulation. The first important field is the attribute ``Cosmological
run`` which is set to ``1`` for cosmological runs and to ``0`` otherwise. This
allows users to quickly distinguish between these two main modes. Most values in
this section only make sense for cosmological runs.

All quantities are expressed in the internal system of units (note that this may
differ from the units used in the particle arrays). Values like the look-back
time are given for the redshift (or scale-factor) of this snapshot.

Policy
~~~~~~

The group ``/Policy`` list the engine policies (defined in ``src/engine.h``)
that were activated in the run that dumped this snapshot. The policies roughly
translate to the main run-time parameters of SWIFT.

GravityScheme
~~~~~~~~~~~~~

HydroScheme
~~~~~~~~~~~

StarsScheme
~~~~~~~~~~~

SubgridScheme
~~~~~~~~~~~~~

This section of the meta-data mainly contains basic information about the
flavour of sub-grid schemes used in the simulation. This is typically a list of
attributes describing the parameters of the model. Users willing to add
information can edit the functions ``chemistry_write_flavour()``,
``cooling_write_flavour()``, etc. located in the i/o header of each scheme.

The other important output stored in that group is the ``NamedColumns``
sub-group. In it, we store the names of the columns of larger particle arrays
that are stored as large n-dimensional arrays. For instance, in the EAGLE model,
the individual chemical element fractions of each particles are stored as a Nx9
array, where N is the number of particles (See
:ref:`EAGLE_chemical_tracers`). This array is labeled ``ElementMassFractions``
and is used instead of 9 individual 1-d arrays. In the ``NamedColumns``
sub-group we store as an array of strings the name of each of the 9 individual
columns. In this case, the name of the 9 elements traced by the model. This
array has the same name as the particle array it corresponds to; here
``ElementMassFractions``. The same mechanism is used for other quantities stored
in a similar fashion. This allows external tools reading SWIFT snapshots to give
meaningful names to more complex entries of the particle arrays.

Unit systems
------------

The snapshots contain *two* groups at the root containing information about the
unit systems used in the snapshots.

The main one ``Units`` contains the units used in the snapshot. In a similar
fashion to what is done for the parameter files (see :ref:`Parameters_units`),
SWIFT specifies only the basic units. These are the units of mass (``U_M``),
length (``U_L``), time (``U_t``), electric current (``U_I``) and temperature
(``U_T``). These are specified in units of their CGS equivalents (gram,
centimeter, second, Ampere, Kelvin). All the quantities present in the particle
arrays are expressed in this system of units. For each quantity, SWIFT gives the
conversion factor in terms of these units. For instance, the internal energy per
unit mass would be expressed as ``U_L^2 U_t^-2``, which in the CGS unit system
translates to :math:`cm/s^2 = erg/g`.

The second group ``InternalCodeUnits`` contains the unit system that was used
internally by the code when running the simulation. This is in most cases the
same system as given in ``Units`` but since users can specify a different
system for the snapshots, there might be cases where they differ. As this system
only relates to what was used inside the code and not in the snapshots
themselves, this group is mostly here to report on the code's run-time behaviour
and is used to express all the quantities in the meta-data (e.g. in the
cosmology group or the softening lengths in the gravity group).

Used and unused run-time parameters
-----------------------------------

The groups ``/Parameters`` and ``UnusedParameters`` located at the root of the file
contain the list of all the run-time parameters used by the run with their
values and the list of parameters that were in the YAML but were not read. The
content of these two groups is identical to the ``used_parameters.yml`` and
``unused_parameters.yml`` files produced by SWIFT when starting a run (See
the :ref:`Parameters_basics` section of the documentation).

Structure of the particle arrays
--------------------------------

There are several groups that contain 'auxiliary' information, such as
``Header``.  Particle data is placed in separate groups depending of the type of
the particles. There are currently 6 particle types available. The type use the
naming convention of Gadget-2 (with the OWLS and EAGLE extensions). A more
intuitive naming convention is given in the form of aliases within the file. The
aliases are shown in the third column of the table.

+---------------------+------------------------+-----------------------------+----------------------------------------+
| HDF5 Group Name     | Physical Particle Type | HDF5 alias                  | In code ``enum part_type``             |
+=====================+========================+=============================+========================================+
| ``/PartType0/``     | Gas                    | ``/GasParticles/``          | ``swift_type_gas``                     |
+---------------------+------------------------+-----------------------------+----------------------------------------+
| ``/PartType1/``     | Dark Matter            | ``/DMParticles/``           | ``swift_type_dark_matter``             |
+---------------------+------------------------+-----------------------------+----------------------------------------+
| ``/PartType2/``     | Background Dark Matter | ``/DMBackgroundParticles/`` | ``swift_type_dark_matter_background``  |
+---------------------+------------------------+-----------------------------+----------------------------------------+
| ``/PartType3/``     | Sinks                  | ``/SinkParticles/``         | ``swift_type_sink``                    |
+---------------------+------------------------+-----------------------------+----------------------------------------+
| ``/PartType4/``     | Stars                  | ``/StarsParticles/``        | ``swift_type_star``                    |
+---------------------+------------------------+-----------------------------+----------------------------------------+
| ``/PartType5/``     | Black Holes            | ``/BHParticles/``           | ``swift_type_black_hole``              |
+---------------------+------------------------+-----------------------------+----------------------------------------+
| ``/PartType6/``     | Neutrino Dark Matter   | ``/NeutrinoParticles/``     | ``swift_type_neutrino``                |
+---------------------+------------------------+-----------------------------+----------------------------------------+

The last column in the table gives the ``enum`` value from ``part_type.h``
corresponding to a given entry in the files.

For completeness, the list of particle type names is stored in the snapshot
header in the array ``/Header/PartTypeNames``. The number of types (aka. the
length of this array) is stored as the attribute ``/Header/NumPartTypes``.

Each group contains a series of arrays corresponding to each field of the
particles stored in the snapshots. The exact list of fields depends on what
compile time options were used and what module was activated. A full list can be
obtained by running SWIFT with the ``-o`` runtime option (See
:ref:`Output_selection_label` for details). Each field contains a short
description attribute giving a brief summary of what the quantity represents.

All the individual arrays created by SWIFT have had the Fletcher 32 check-sum
filter applied by the HDF5 library when writing them. This means that any
eventual data corruption on the disks will be detected and reported by the
library when attempting to read the data.

Additionally, some compression filter may have been applied to the fields. See
the :ref:`Parameters_snapshots` section of the parameter file description for
more details.

Unit information for individual fields
--------------------------------------

Each particle field contains meta-data about the units and how to
convert it to CGS in physical or co-moving frames. The meta-data is in
part designed for users to directly read and in part for machine
reading of the information. Each field contains the exponent of the
scale-factor, reduced Hubble constant [#f2]_ and each of the 5 base units
that is required to convert the field values to physical CGS
units. These fields are:

+----------------------+---------------------------------------+
| Meta-data field name | Description                           |
+======================+=======================================+
| ``U_L exponent``     | Power of the length unit              |
+----------------------+---------------------------------------+
| ``U_M exponent``     | Power of the mass unit                |
+----------------------+---------------------------------------+
| ``U_t exponent``     | Power of the time unit                |
+----------------------+---------------------------------------+
| ``U_I exponent``     | Power of the current unit             |
+----------------------+---------------------------------------+
| ``U_T exponent``     | Power of the temperature unit         |
+----------------------+---------------------------------------+
| ``a-scale exponent`` | Power of the scale-factor             |
+----------------------+---------------------------------------+
| ``h-scale exponent`` | Power of the reduced Hubble constant  |
+----------------------+---------------------------------------+

These are used by the ``swiftsimio`` python library to read units and
we encourage users to use this meta-data directly in their automated
tools.

As an example, the fluid densities (which are written in the co-moving
frame) have the following (non-zero) conversion factors:

 * ``U_L exponent``: -3
 * ``U_M exponent``: 1
 * ``a-scale exponent``: -3

This condensed information is stored in the string ``Expression for
physical CGS units``, which in the case of the densities would read
``a^-3 U_M U_L^-3 [ g cm^-3 ]``. The values of the ``U_x`` can be
found in the ``Units System`` group at the root of the snapshot (see
above). Note that only unit factors with non-zero exponents are
printed to this string.

Additionally, the meta-data contains the numerical conversion factor
from the field to co-moving CGS and physical CGS assuming the units in
the ``Unit System`` group. These are:

 * ``Conversion factor to CGS (not including cosmological corrections``
 * ``Conversion factor to phyical CGS (including cosmological corrections)``

These are designed for the users to directly use if they don't want to
compute the individual exponents themselves. As an example, in the
case of the densities and assuming the usual system of units
(:math:`10^{10} \rm{M}_\odot`, :math:`100 \rm{km/s}`, :math:`\rm{Mpc}`) at redshift
0.1, the conversion factors are:

 * Conversion to CGS: :math:`6.76814403 \times 10^{-31}`
 * Conversion to physical CGS: :math:`9.00808555 \times 10^{-31}`

In the case of a non-cosmological simulation, these two expressions
are identical since :math:`a=1`.

Particle splitting metadata
---------------------------

When particle splitting is turned on (see :ref:`Parameters_basics`; by using
``particle_splitting=1`` in the parameter file) some particles in the output
may have been created from the 'splitting' of a single, over-massive, particle.

There are three fields, associated with all gas, star, and black hole particles,
that can be used to understand if, and how, these particles were split.

These three fields are:

+ ``ProgenitorIDs``, the IDs of the gas particles in the initial conditions
  that is the direct progenitor of this particle.
+ ``SplitCounts``, the number of times this gas particle has been split; or,
  if a star or black hole, how many times the gas particle that became this
  star (or black hole seed) was split before becoming so.
+ ``SplitTrees``, a binary tree (encoded as a 64 bit integer) showing how this
  particle was split. Each item in the tree shows whether this particle retained
  its original ID (encoded as 0) or was given a new ID (encoded as 1) in the
  splitting event. This data is enough to completely reconstruct the splitting 
  history of the particles.

For example, if a particle has been split 5 times (``SplitCounts=5`` for this
particle), and has a binary tree of "10010", it retained its original ID in
the first event, was given a new one in the second event, for the next two
events it retained its new ID (obtained in the second event), and finally was
given a new ID in the final event. Throughout this process, the value of
``ProgenitorIDs`` remained the same. Through this system, we can ensure that
the combination of ``ProgenitorID`` and this binary tree corresponds to a
fully traceable, unique, identifier for every particle in the simulation volume.

Note that we can only track 64 splitting events for a given particle, and after
this the binary tree is meaningless. In practice, however, such a high number
of splitting events is extremely unlikely to occur.

An example is provided in ``examples/SubgridTests/ParticleSplitting``, with
a figure showing how one particle is split (eventually) into 16 descendants
that makes use of this metadata.
   
Quick access to particles via hash-tables
-----------------------------------------

The particles are not sorted in a specific order when they are written to the
snapshots. However, the particles are sorted into the top-level cell structure
used internally by the code every time a tree rebuild is triggered. The
top-level cells are a coarse-grained mesh but knowing which particle belongs to
which cell can nevertheless be useful to rapidly access particles in a given
region only.

One important caveat is that particles are free to drift out of their cells
between rebuilds of the tree (but not by more than one cell-length). If one
wants to have all the particles in a given cell, one has to read all the
neighbouring cells as well. We note that for image making purposes, for instance
to generate a slice, this is typically not necessary and reading just the cells
of interest is sufficient.

At the root of the HDF5 file, the ``Cells`` group contains all the relevant
information. The dimension of the top-level grid (a triplet of integers) is
given by the attribute ``Cells/Meta-data/dimension`` and the size of each cell (a
triplet of floating-point numbers) is given by the attribute
``Cells/Meta-data/size``. All the cells have the same size but for non-cubic
simulation volumes the cells themselves can have different sizes along each
axis.

The ``/Cells/Centres`` array gives the centre of each of the top-level cells in
the simulation volume. Both the cell sizes and positions of the centres are
expressed in the unit system used for the snapshots (see above) and are hence
consistent with the particle positions themselves. 

Once the cell(s) containing the region of interest has been located,
users can use the ``/Cells/Files/PartTypeN/``,
``/Cells/Counts/PartTypeN/`` and
``/Cells/OffsetsInFile/PartTypeN/`` to retrieve the location of
the particles of type ``N`` in the ``/PartTypeN`` arrays.  These
contain information about which file contains the particles of a given
cell. It also gives the offset from the start of the ``/PartTypeN``
array *in that file* at which the particles of that cell are located
and how many particles are in the cell. This allows to read a single
contiguous section of the whole array by directly reading the slab
starting at the offset and with the given length.

The cells, files, offsets in file and counts arrays are sorted
spatially using C-style ordering. That means the inner-most loop runs
over the z axis, then y axis and x is the slowest varying dimension.

In the case of a single-file snapshot, the ``Files`` array is just an array of
zeroes since all the particles will be in the 0-th file. Note also that in the
case of a multi-files snapshot, a cell is always contained in a single file.

As noted above, particles can (slightly) drift out of their cells. This can be
problematic in cases where one wants to find precisely all the particles in a
given region. To help with this, the meta-data also contains a "cell bounding
box". The arrays ``/Cells/MinPositions/PartTypeN`` and
``/Cells/MaxPositions/PartTypeN`` contain the minimal (maximal) x,y,z
coordinates of all the particles of this type in the cells. Note that these
coordinates can be outside of the cell itself. When using periodic boundary
conditions, no box-wrapping is applied.

If a snapshot used a sub-sampled output, then the counts and offsets are
adjusted accordingly and correspond to the actual content of the file
(i.e. after the sub-sampling was applied).

As an example, if one is interested in retriving all the densities of the gas
particles in the cell around the position `[1, 1, 1]` in a single-file
snapstshot one could use a piece of code similar to:

.. code-block:: python
   :linenos:

   import numpy as np
   import h5py

   snapshot_file = h5py.File("snapshot.hdf5", "r")

   my_pos = [1, 1, 1]

   # Read in the cell centres and size
   nr_cells = f["/Cells/Meta-data"].attrs["nr_cells"]
   centres = f["/Cells/Centres"][:,:]
   size = f["/Cells/Meta-data"].attrs["size"]
   half_size = size / 2.

   # Look for the cell containing the position of interest.
   #
   # Note that since the cells are sorted spatially, we would formally
   # not need to do this search and could jump directly to the correct 'i'.
   my_cell = -1
   for i in range(nr_cells):
      if my_pos[0] > centres[i, 0] - half_size[0] and my_pos[0] < centres[i, 0] + half_size[0] and
         my_pos[1] > centres[i, 1] - half_size[1] and my_pos[1] < centres[i, 1] + half_size[1] and
         my_pos[2] > centres[i, 2] - half_size[2] and my_pos[2] < centres[i, 2] + half_size[2]:
	 my_cell = i
	 break
   
   # Print the position of the centre of the cell of interest
   centre = snapshot_file["/Cells/Centres"][my_cell, :]
   print("Centre of the cell:", centre)

   # Retrieve the offset and counts
   my_offset = snapshot_file["/Cells/OffsetsInFile/PartType0"][my_cell]
   my_count = snapshot_file["/Cells/Counts/PartType0"][my_cell]

   # Get the densities of the particles in this cell
   rho = snapshot_file["/PartType0/Density"][my_offset:my_offset + my_count]

For large simulations, this vastly reduces the amount of data that needs to be read
from the disk.

Note that this is all automated in the ``swiftsimio`` python library
and we highly encourage its use.

Meta-file for distributed snapshots
-----------------------------------

If distributed snapshots are chosen for an MPI parallel run (see
:ref:`Parameters_snapshots`), N snapshot files are produced, where N is the
number of MPI ranks. When HDF5 1.10.0 or higher is available, an
additional meta-snapshot is produced that uses HDF5's virtual dataset
feature to present these N files as if they were a single, regular
snapshot file.

The meta-snapshot contains all the meta-data (including the top level
cell hash-tables) contained in a regular snapshot, but does not store
any actual particle data. Instead, the particle datasets contain virtual
links to the corresponding particle data in the distributed snapshot
files. Since this is a feature of the HDF5 library itself, this is
entirely transparent to modules like ``h5py`` that try to read the data.
A user only needs to access the meta-snapshot, and the HDF5 library
takes care of the rest.

The virtual links in the meta-snapshot only work if the HDF5 library
knows the location of the distributed snapshots. These are stored within
the meta-snapshot as relative paths. When SWIFT produces a distributed
snapshot, all files are placed within the same directory. This means
that the meta-snapshot can only be safely read if the other N files are
also present in the same directory.

The header of a meta-snapshot looks exactly like the header of a normal,
non-distributed snapshot (i.e. ``NumFilesPerSnapshot`` is 1). However,
the attribute ``Virtual`` is set to 1 to distinguish it from a normal
snapshot file.

.. [#f1] In the rare case where an output
	 selection (see :ref:`Output_selection_label`) disabling a given particle type in
	 its entirety was used, the corresponding entry in ``NumPart_ThisFile`` will be 0
	 whilst the ``NumPart_Total`` field will still contain the number of
	 particles present in the run.


.. [#f2] Note that all quantities in SWIFT are always "h-free" in the sense that
	 they are expressed in units withouy any h terms. This implies that the
	 ``h-scale exponent`` field value is always 0. SWIFT nevertheless
	 includes this field to be comprehensive and to prevent confusion with
	 other software packages that express their quantities with h-full
	 units.
