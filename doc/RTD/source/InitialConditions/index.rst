.. Initial Conditions
   Josh Borrow, 5th April 2018

Initial Conditions
==================

To run anything more than examples from our suite, you will need to be able to
produce your own initial conditions for SWIFT. We use the same initial
conditions format as the popular `GADGET-2
<https://wwwmpa.mpa-garching.mpg.de/~volker/gadget/>`_ code, which uses HDF5 for
its type 3 format. Note that we do not support the GADGET-2 types 1 and 2
formats.

One crucial difference is that whilst GADGET-2 can have initial conditions split
over many files SWIFT only supports initial conditions in one single file. **ICs
split over multiple files cannot be read by SWIFT**. See the
":ref:`multiple_files_ICs`" section below for possible solutions. In GADGET-2
having multiple files allows multiple ones to be read in parallel and is the
only way the code can handle more than 2^31 particles. This limitation is not in
place in SWIFT. A single file can contain any number of particles (well... up to
2^64...)  and the file is read in parallel by HDF5 when running on more than one
compute node.

The original GADGET-2 file format only contains 2 types of particles: gas
particles and 5 sorts of collision-less particles that allow users to run with 5
separate particle masses and softenings. In SWIFT, we expand on this by using
two of these types for stars and black holes.

As the original documentation for the GADGET-2 initial conditions format is
quite sparse, we lay out here all of the necessary components. If you are
generating your initial conditions from python, we recommend you use the h5py
package. We provide a writing wrapper for this for our initial conditions in
``examples/KeplerianRing/write_gadget.py``.

You can find out more about the HDF5 format on their `webpages
<https://support.hdfgroup.org/HDF5/doc/H5.intro.html>`_.


Structure of the File
---------------------

There are several groups that contain 'auxiliary' information, such as
``Header``.  Particle data is placed in separate groups depending of the type of
the particles. Some types are currently ignored by SWIFT but are kept in the
file format for compatibility reasons.

+---------------------+------------------------+----------------------------+
| HDF5 Group Name     | Physical Particle Type | In code ``enum part_type`` |
+=====================+========================+============================+
| ``/PartType0/``     | Gas                    | ``swift_type_gas``         |
+---------------------+------------------------+----------------------------+
| ``/PartType1/``     | Dark Matter            | ``swift_type_dark_matter`` |
+---------------------+------------------------+----------------------------+
| ``/PartType2/``     | Ignored                |                            |
+---------------------+------------------------+----------------------------+
| ``/PartType3/``     | Ignored                |                            |
+---------------------+------------------------+----------------------------+
| ``/PartType4/``     | Stars                  | ``swift_type_star``        |
+---------------------+------------------------+----------------------------+
| ``/PartType5/``     | Black Holes            | ``swift_type_black_hole``  |
+---------------------+------------------------+----------------------------+

The last column in the table gives the ``enum`` value from ``part_type.h``
corresponding to a given entry in the files.

Note that the only particles that have hydrodynamical forces calculated between
them are those in ``PartType0``.


Necessary Components
--------------------

There are several necessary components (in particular header information) in a
SWIFT initial conditions file. Again, we recommend that you use the ``write_gadget``
script.

Header
~~~~~~

In the ``/Header/`` group, the following attributes are required:

+ ``Dimension``, an integer indicating the dimensionality of the ICs (1,2 or 3).
  Note that this parameter is an addition to the GADGET-2 format and will be
  ignored by GADGET. SWIFT will use this value to verify that the dimensionality
  of the code matches the ICs. If this parameter is not provided, it defaults
  to 3.
+ ``BoxSize``, a floating point number or N-dimensional (usually 3) array that
  describes the size of the box. If only one number is provided (as per the
  GADGET-2 standard) then the box is assumed have the same size along all the
  axis. In cosmological runs, this is the comoving box-size expressed in the
  units specified in the ``/Units`` group (see below). Note that, unlike GADGET,
  we express all quantities in "h-free" units. So that, for instance, we express
  the box side-length in ``Mpc`` and not ``Mpc/h``. 
+ ``NumPart_Total``, a length 6 array of integers that tells the code how many
  particles of each type are in the initial conditions file. Unlike traditional
  GADGET-2 files, these can be >2^31.
+ ``NumPart_Total_HighWord``, a historical length-6 array that tells the code
  the number of high word particles in the initial conditions there are. If you
  are unsure, just set this to ``[0, 0, 0, 0, 0, 0]``. This does have to be
  present but can be a set of 0s unless you have more than 2^31 particles and
  want to be fully compliant with GADGET-2. Note that, as SWIFT supports
  ``NumPart_Total`` to be >2^31, the use of ``NumPart_Total_HighWord`` is only
  here for compatibility reasons.
+ ``Flag_Entropy_ICs``, a historical value that tells the code if you have
  included entropy or internal energy values in your initial conditions files.
  Acceptable values are 0 or 1. We recommend using internal energies over
  entropy in the ICs and hence have this flag set to 0.

You may want to include the following for backwards-compatibility with many
GADGET-2 based analysis programs:

+ ``MassTable``, an array of length 6 which gives the masses of each particle
  type. SWIFT ignores this and uses the individual particle masses, but some
  programs will crash if it is not included.
+ ``NumPart_ThisFile``, a length 6 array of integers describing the number of
  particles in this file. If you have followed the above advice, this will be
  exactly the same as the ``NumPart_Total`` array. As SWIFT only uses ICs
  contained in a single file, this is not necessary for SWIFT-only ICs.
+ ``NumFilesPerSnapshot``, again a historical integer value that tells the code
  how many files there are per snapshot. You will probably want to set
  this to 1. If this field is present in a SWIFT IC file and has a
  value different from 1, the code will return an error message.
+ ``Time``, time of the start of the simulation in internal units or expressed
  as a scale-factor for cosmological runs. SWIFT ignores this and reads it from
  the parameter file.


Particle Data
~~~~~~~~~~~~~

Now for the interesting part! You can include particle data groups for each
individual particle type (e.g. ``/PartType0/``) that have the following *datasets*:

+ ``Coordinates``, an array of shape (N, 3) where N is the number of particles
  of that type, that are the cartesian co-ordinates of the
  particles. Co-ordinates must be within the box so, in the case of a cube
  within [0, L)^3 where L is the side-length of the simulation volume. In the
  case of cosmological simulations, these are the co-moving positions.
+ ``Velocities``, an array of shape (N, 3) that is the cartesian velocities of
  the particles. When running cosmological simulations, these are the peculiar
  velocities. Note that this is different from GADGET which uses peculiar
  velocities divided by ``sqrt(a)`` (see below for a fix).
+ ``ParticleIDs``, an array of length N that are unique identifying numbers for
  each particle. Note that these have to be unique to a particle, and cannot be
  the same even between particle types. The **IDs must be >= 0**. Negative
  IDs will be rejected by the code.
+ ``Masses``, an array of length N that gives the masses of the particles.

For ``PartType0`` (i.e. particles that interact through hydro-dynamics), you will
need the following auxiliary items:

+ ``SmoothingLength``, the smoothing lengths of the particles. These will be
  tidied up a bit, but it is best if you provide accurate numbers. In
  cosmological runs, these are the co-moving smoothing lengths.
+ ``InternalEnergy``, an array of length N that gives the internal energies per
  unit mass of the particles. If the hydro-scheme used in the code is based on
  another thermodynamical quantity (entropy or total energy, etc.), the
  conversion will happen inside the code. In cosmological runs, this is the
  **physical** internal energy per unit mass. This has the dimension of velocity
  squared.

  
Note that for cosmological runs, all quantities have to be expressed in "h-free"
dimensions. This means ``Mpc`` and not ``Mpc/h`` for instance. If the ICs have
been generated for GADGET (where h-full values are expected), the parameter
``InitialConditions:cleanup_h_factors`` can be set to ``1`` in the
:ref:`Parameter_File_label` to make SWIFT convert the quantities read in to
h-free quantities. Switching this parameter on will also affect the box size
read from the ``/Header/`` group (see above).

Similarly, GADGET cosmological ICs have traditionally used velocities expressed
as peculiar velocities divided by ``sqrt(a)``. This can be undone by switching
on the parameter ``InitialConditions:cleanup_velocity_factors`` in the
:ref:`Parameter_File_label`.


.. _ICs_units_label:

Optional Components
-------------------

In the ``/Units/`` HDF5 group, you cans specify what units your initial conditions are
in. If this group is not present, the code assumes that you are using the same
units for your initial conditions as in your :ref:`Parameter_File_label`
(i.e. as the internal units system used by the code), but it is best to include
them to be on the safe side. You will need:

+ ``Unit length in cgs (U_L)``
+ ``Unit mass in cgs (U_M)``
+ ``Unit time in cgs (U_t)``
+ ``Unit current in cgs (U_I)``
+ ``Unit temperature in cgs (U_T)``

These are all floating point numbers. Note that we specify the time units and
not the velocity units.

If the units specified in the initial conditions are different from the internal
units (specified in the parameter file), SWIFT will perform a conversion of all
the quantities when reading in the ICs. This includes a conversion of the box
size read from the ``/Header/`` group.


     
Summary
-------

You should have an HDF5 file with the following structure:

.. code-block:: bash

   Header/
     BoxSize=[x, y, z]
     Flag_Entropy_ICs=0
     NumPart_Total=[0, 1, 0, 0, 4, 5]
     NumPart_Total_HighWord=[0, 0, 0, 0, 0, 0]
   Units/
     Unit current in cgs (U_I)=1.0
     Unit length in cgs (U_L)=1.0
     Unit mass in cgs (U_M)=1.0
     Unit temperature in cgs (U_T)=1.0
     Unit time in cgs (U_t)=1.0
   PartType0/
     Coordinates=[[x, y, z]]
     Velocities=[[vx, vy, vz]]
     ParticleIDs=[...]
     Masses=[...]
     InternalEnergy=[...]
     SmoothingLength=[...]
   PartType1/
     Coordinates=[[x, y, z]]
     Velocities=[[vx, vy, vz]]
     ParticleIDs=[...]
     Masses=[...]

.. _multiple_files_ICs:
     
ICs split over multiple files
-----------------------------

A basic script ``tools/combine_ics.py`` is provided to merge basic GADGET-2
initial conditions split into multiple files into one single valid file. This
script can handle simple HDF5 files (GADGET-2 type 3 ICs) that follow the format
described above but split over multiple files.

The script can also convert ICs using a ``MassTable`` and create the
corresponding particle fields. Note that additional fields present in ICs beyond
the simple GADGET-2 specification will not be merged.

One additional option is to compress the fields in the files using HDF5's gzip
compression. This is very effective for the fields such as masses or particle
IDs which are very similar. A checksum filter is also applied in all cases to
help with data curation.

**We caution that this script is very basic and should only be used with great
caution.** 



