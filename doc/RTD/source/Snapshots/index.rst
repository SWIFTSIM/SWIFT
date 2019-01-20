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
the particles. The type use the naming convention of Gadget-2 (with
the OWLS and EAGLE extensions).

+---------------------+------------------------+----------------------------+
| HDF5 Group Name     | Physical Particle Type | In code ``enum part_type`` |
+=====================+========================+============================+
| ``/PartType0/``     | Gas                    | ``swift_type_gas``         |
+---------------------+------------------------+----------------------------+
| ``/PartType1/``     | Dark Matter            | ``swift_type_dark_matter`` |
+---------------------+------------------------+----------------------------+
| ``/PartType4/``     | Stars                  | ``swift_type_star``        |
+---------------------+------------------------+----------------------------+
| ``/PartType5/``     | Black Holes            | ``swift_type_black_hole``  |
+---------------------+------------------------+----------------------------+

The last column in the table gives the ``enum`` value from ``part_type.h``
corresponding to a given entry in the files.

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

The ``/Cells/Centres`` array gives the centre of each of the top-level cells in the
simulation volume. Both the cell sizes and positions of the centres are
expressed in the unit system used for the snapshots (see above) and are hence
consistent with the particle positions themselves.

Once the cell(s) containing the region of interest has been located, users can
use the ``/Cells/Offsets/PartTypeN/Counts`` and
``/Cells/Offsets/PartTypeN/Offsets`` to retrieve the location of the particles
of type ``N`` in the ``/PartTypeN`` arrays. For instance, if one is interested
in retriving all the densities of the gas particles in the cell around the
position `[1, 1, 1]` one could use a piece of code similar to:

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

   # Look for the cell containing the position of interest
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
   my_offset = snapshot_file["/Cells/Offsets/PartType0"][my_cell]
   my_count = snapshot_file["/Cells/Counts/PartType0"][my_cell]

   # Get the densities of the particles in this cell
   rho = snapshot_file["/PartType0/Density"][my_offset:my_offset + my_count]

For large simulations, this vastly reduces the amount of data that needs to be read
from the disk.
