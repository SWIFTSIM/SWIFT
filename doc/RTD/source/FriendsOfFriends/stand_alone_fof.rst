.. Friends Of Friends
   Matthieu Schaller 15th June 2019

.. _fof_stand_alone_label:

Stand-Alone Friends-Of-Friends
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The stand-alone version of the Friends-Of-Friends module can be
compiled by configuring the code with the option
``--enable-stand-alone-fof``. The ``fof`` and ``fof_mpi`` executables
will then be generated alongside the regular SWIFT ones.

The executable takes a parameter file as an argument. It will then read the
snapshot specified in the parameter file (specified as an initial condition
file) and extract all the dark matter particles by default. FOF is then run on
these particles and a catalogue of groups is written to disk. Additional
particle types can be read and processed by the stand-alone FOF code by adding
any of the following runtime parameters to the command line:

 * ``--hydro``: Read and process the gas particles,
 * ``--stars``: Read and process the star particles,
 * ``--black-holes``: Read and process the black hole particles,
 * ``--sinks``: Read and process the sink particles,
 * ``--cosmology``: Consider cosmological terms.

Running with cosmology is necessary when using a linking length based
on the mean inter-particle separation. For analysis where an absolute
linking length is used this is not necessary.

The group catalogues contain a unique identifier for each group as
well as the number of particles in each group and their total mass (in
internal units). The FOF code will also write a snapshot with an
additional field for each particle. This contains the ``GroupID`` of
each particle and can be used to find all the particles in a given
halo and to link them to the information stored in the catalogue.

.. warning::

   Note that since not all particle properties are read in stand-alone
   mode, not all particle properties will be written to the snapshot generated
   by the stand-alone FOF.

