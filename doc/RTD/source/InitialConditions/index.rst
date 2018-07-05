.. Initial Conditions
   Josh Borrow, 5th April 2018

Initial Conditions
==================

To run anything more than examples from our suite, you will need to be able to 
produce your own initial conditions for SWIFT. We use the same initial conditions
as the popular GADGET-2 code, which uses the HDF5 file format.

As the original documentation for the GADGET-2 initial conditions format is
quite sparse, we lay out here all of the necessary components. If you are generating
your initial conditions from python, we recommend you use the h5py package. We
provide a writing wrapper for this for our initial conditions in
``examples/KeplerianRing/write_gadget.py``.

You can find out more about the HDF5 format on their webpages, here:
https://support.hdfgroup.org/HDF5/doc/H5.intro.html


Structure of the File
---------------------

There are several groups that contain 'auxilliary' information, such as ``Header``.
Particle data is placed in groups that signify particle type.

+---------------------+------------------------+
| Group Name          | Physical Particle Type |
+=====================+========================+
| ``PartType0``       | Gas                    |
+---------------------+------------------------+
| ``PartType1``       | Dark Matter            |
+---------------------+------------------------+
| ``PartType2``       | Ignored                |
+---------------------+------------------------+
| ``PartType3``       | Ignored                |
+---------------------+------------------------+
| ``PartType4``       | Stars                  |
+---------------------+------------------------+
| ``PartType5``       | Black Holes            |
+---------------------+------------------------+

Currently, not all of these particle types are included in SWIFT. Note that the
only particles that have hydrodynamical forces calculated between them are those
in ``PartType0``.


Necessary Components
--------------------

There are several necessary components (in particular header information) in a
SWIFT initial conditions file. Again, we recommend that you use the ``write_gadget``
script.

Header
~~~~~~

In ``Header``, the following attributes are required:

+ ``BoxSize``, a floating point number or N-dimensional (usually 3) array
  that describes the size of the box.
+ ``Flag_Entropy_ICs``, a historical value that tells the code if you have
  included entropy or internal energy values in your intial conditions files.
  Acceptable values are 0 or 1.
+ ``NumPart_Total``, a length 6 array of integers that tells the code how many
  particles are of each type are in the initial conditions file.
+ ``NumPart_Total_HighWord``, a historical length-6 array that tells the code 
  the number of high word particles in the initial conditions there are. If
  you are unsure, just set this to ``[0, 0, 0, 0, 0, 0]``. This does have to be
  present, but unlike GADGET-2, this can be a set of 0s unless you have more than
  2^31 particles.
+ ``NumFilesPerSnapshot``, again a historical integer value that tells the code
  how many files there are per snapshot. You will probably want to set this to 1
  and simply have a single HDF5 file for your initial conditions; SWIFT can
  leverage parallel-HDF5 to read from this single file in parallel.
+ ``NumPart_ThisFile``, a length 6 array of integers describing the number of
  particles in this file. If you have followed the above advice, this will be
  exactly the same as the ``NumPart_Total`` array.

You may want to include the following for backwards-compatibility with many
GADGET-2 based analysis programs:

+ ``MassTable``, an array of length 6 which gives the masses of each particle
  type. SWIFT ignores this and uses the individual particle masses, but some
  programs will crash if it is not included.
+ ``Time``, the internal code time of the start (set this to 0).

RuntimePars
~~~~~~~~~~~

In ``RuntimePars``, the following attributes are required:

+ ``PeriodicBoundaryConditionsOn``, a flag to tell the code whether or not you
  have periodic boundaries switched on. Again, this is historical; it should be
  set to 1 (default) if you have the code running in periodic mode, or 0 otherwise.


Units
~~~~~

In ``Units``, you will need to specify what units your initial conditions are
in. If these are not present, the code assumes that you are using the same
units for your initial conditions as are in your parameterfile, but it is best
to include them to be on the safe side. You will need:

+ ``Unit current in cgs (U_I)``
+ ``Unit length in cgs (U_L)``
+ ``Unit mass in cgs (U_M)``
+ ``Unit temperature in cgs (U_T)``
+ ``Unit time in cgs (U_t)``

These are all floating point numbers.


Particle Data
~~~~~~~~~~~~~

Now for the interesting part! You can include particle data groups for each
individual particle type (e.g. ``PartType0``) that have the following _datasets_:

+ ``Coordinates``, an array of shape (N, 3) where N is the number of particles
  of that type, that are the cartesian co-ordinates of the particles. Co-ordinates
  must be positive, but will be wrapped on reading to be within the periodic box.
+ ``Velocities``, an array of shape (N, 3) that is the cartesian velocities 
  of the particles.
+ ``ParticleIDs``, an array of length N that are unique identifying numbers for
  each particle. Note that these have to be unique to a particle, and cannot be
  the same even between particle types. Please ensure that your IDs are positive
  integer numbers.
+ ``Masses``, an array of length N that gives the masses of the particles.

For ``PartType0`` (i.e. particles that interact through hydrodynamics), you will
need the following auxilliary items:

+ ``InternalEnergy``, an array of length N that gives the internal energies of
  the particles. For PressureEntropy, you can specify ``Entropy`` instead.
+ ``SmoothingLength``, the smoothing lenghts of the particles. These will be
  tidied up a bit, but it is best if you provide accurate numbers.


Summary
~~~~~~~

You should have an HDF5 file with the following structure:

.. code-block:: bash

   Header/
     BoxSize=[x, y, z]
     Flag_Entropy_ICs=1
     NumPart_Total=[0, 1, 2, 3, 4, 5]
     NumPart_Total_HighWord=[0, 0, 0, 0, 0, 0]
     NumFilesPerSnapshot=1
     NumPart_ThisFile=[0, 1, 2, 3, 4, 5]
   RuntimePars/
     PeriodicBoundariesOn=1
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


