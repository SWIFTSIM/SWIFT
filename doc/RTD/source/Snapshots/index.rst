.. Snapshots
   Matthieu Schaller, 5th January 2019

.. _snapshots:

Snapshots
=========

The snapshots are stored using the HDF5 format and are fully compatible with
Gadget-2. They do, however, contain a large set of extensions including units,
meta-data about the code and runs as well as facilities to quickly access the
particles in a specific region of the simulation volume.

Meta-data about the code and run
--------------------------------

Unit system
-----------

Used and unused run-time parameters
-----------------------------------

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
