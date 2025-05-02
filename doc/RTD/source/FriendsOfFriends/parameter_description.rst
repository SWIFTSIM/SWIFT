.. Friends Of Friends
   Matthieu Schaller 15th June 2019

.. _Fof_Parameter_Description_label:

Friends-Of-Friends Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``FOF`` section of the parameter file contains all the options related
to the group finder. Some parameters only make sense in the stand-alone
version and some only in the on-the-fly version of the module.

The main parameter is the linking length that will be used. The easiest way
to define it is to set the ratio of the linking length to the mean *dark
matter* inter-particle separation (irrespective of whether there are other
kinds of particles present in the run). This is done using the parameter
``linking_length_ratio`` and the typical value used is
``0.2``. Users can optionally overwrite the linking length by imposing an
absolute value using the parameter ``absolute_linking_length``. This is
expressed in internal units. This value will be ignored (and the ratio of
the mean inter-particle separation will be used) when set to ``-1``.

The categories of particles are specified using the ``linking_types`` and
``attaching_types`` arrays. They are of the length of the number of particle
types in SWIFT (currently 7) and specify for each type using ``1`` or ``0``
whether or not the given particle type is in this category. Types not present
in either category are ignored entirely.

The second important parameter is the minimal size of groups to retain in
the catalogues. This is given in terms of number of particles (of all types)
via the parameter ``min_group_size``. When analysing simulations, to
identify haloes, the common practice is to set this to ``32`` in order to
not plague the catalogue with too many small, likely unbound, structures.
When using the FOF code on-the-fly for black hole seeding, larger values
are recommended as there is no need to store groups much smaller than the
minimal halo mass (see below).

------------------------

In the case of black holes seeding, we run the FOF module on-the-fly during
a cosmological simulation. Black hole seeding is enabled with the following
parameter:

  * Enable seeding of black holes in FOF groups: ``seed_black_holes_enabled``

This should be set to 1 to enable or 0 to disable black hole seeding. The
time of the first FOF call for seeding is controlled by the following two
options:

  * Time of the first FOF call (non-cosmological runs): ``time_first``,
  * Scale-factor of the first FOF call (cosmological runs):
    ``scale_factor_first``.

One of those two parameters has to be provided depending on the type of
run. In the case of non-cosmological runs, the time of the first FOF call
is expressed in the internal units of time. Users also have to provide the
difference in time (or scale-factor) between consecutive outputs:

  * Time difference between consecutive outputs: ``delta_time``.

The last parameter relevant to the on-the-fly FOF module is the minimal
halo (group) mass considered for seeding black holes. This is specified by
the parameter ``black_hole_seed_halo_mass_Msun`` which is expressed in
solar masses.

There are two ways to invoke FOF on the fly for purposes other than black hole
seeding:

Firstly, one can switch on the ``invoke_fof`` parameter in the
``Snapshots`` section of the parameter file. This will produce a catalogue every
time the code writes a snapshot.

The second option is to set ``dump_catalogue_when_seeding`` in the ``FOF``
section. This will force the code to write a catalogue every time the BH seeding
code is run. 

------------------------

In the case of the stand-alone module, the five seeding parameters
described above are ignored but an additional one needs to be
specified. This is the name of the file in which the catalogue of groups will
be written. This is specified by the parameter ``fof_output``. The linking
length and minimal group to keep are set in the way described above.

There are two additional optional parameters that can be used to modify
details of the stand-alone FOF code. The ``GroupID`` of the particles that
do not belong to any groups is by default set to 2147483647
(i.e. :math:`2^{31}-1`). A different value can be specified by changing the
optional parameter ``group_id_default``. Similarly, the first group in the
catalogue (i.e. the largest group) carries the ``GroupID`` 1. This can be
changed by tweaking the optional parameter ``group_id_offset``.


------------------------

A full FOF section of the YAML parameter file looks like:

.. code:: YAML

    # Parameters of the Friends-Of-Friends module
    FOF:
       basename:                        fof_output  # Filename for the FOF outputs.
       scale_factor_first:              0.91        # Scale-factor of first FoF black hole seeding calls.
       time_first:                      0.2         # Time of first FoF black hole seeding calls.
       delta_time:                      1.005       # Time between consecutive FoF black hole seeding calls.
       min_group_size:                  256         # The minimum no. of particles required for a group.
       linking_types:                   [0, 1, 0, 0, 0, 0, 0]  # Which particle types to consider for linking    (here only DM)
       attaching_types:                 [1, 0, 0, 0, 1, 1, 0]  # Which particle types to consider for attaching  (here gas, stars, and BHs)
       linking_length_ratio:            0.2         # Linking length in units of the main inter-particle separation.
       seed_black_holes_enabled:        0           # Do not seed black holes when running FOF
       black_hole_seed_halo_mass_Msun:  1.5e10      # Minimal halo mass in which to seed a black hole (in solar masses).
       dump_catalogue_when_seeding:     0           # (Optional) Write a FOF catalogue when seeding black holes. Defaults to 0 if unspecified.
       absolute_linking_length:         -1.         # (Optional) Absolute linking length (in internal units).
       group_id_default:                2147483647  # (Optional) Sets the group ID of particles in groups below the minimum size.
       group_id_offset:                 1           # (Optional) Sets the offset of group ID labelling. Defaults to 1 if unspecified.
