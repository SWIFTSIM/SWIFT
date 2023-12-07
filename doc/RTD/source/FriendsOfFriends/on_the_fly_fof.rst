.. Friends Of Friends
   Matthieu Schaller 15th June 2019

.. _fof_on_the_fly_label:

On-the-fly Friends-Of-Friends and Black Hole seeding
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The main purpose of the on-the-fly FOF is to identify haloes during a
cosmological simulation in order to seed some of them with black holes
based on physical considerations.

.. warning::
   In this mode, no group catalogue is written to the disk. The resulting list
   of haloes is only used internally by SWIFT.

Note that a catalogue can nevertheless be written after every seeding call by
setting the optional parameter ``dump_catalogue_when_seeding``.

Once the haloes have been identified by the FOF code, SWIFT will iterate
over the list of groups and will check whether each halo obeys the
following criteria:

  * Is above a user-specified mass threshold (typically
    :math:`10^{10}~\rm{M}_\odot` or more).
  * Contains *at least* one gas particle.
  * Does *not* contain any already existing black hole particle.

If a group satisfies all these requirements then the *densest* gas particle
(based on their SPH density) in the group will be converted to a black
hole. In practice this means that the gas particle is removed and a black
hole particle with the same dynamical mass is inserted at the same position
in the domain and with the same velocity. Additional properties are copied
from the gas particle to the black hole particle depending on the specifics
of the sub-grid model in use (see :ref:`subgrid`).

Given that the halo mass considered for seeding black holes is usually many
hundred times the mass of a single particle, it can be advantageous to
raise the minimal group length required for a halo to appear in the catalogue
of objects. This reduces the number of groups that are kept in memory and
speeds up the seeding procedure by removing haloes that are obviously too
small. For instance, in the case of EAGLE-like runs, the dark matter
particle mass is around :math:`10^7\rm~{M}_\odot` and the minimal halo mass
considered for seeding a black hole is of order
:math:`10^{10}~\rm{M}_\odot`. Groups will hence need to have at least 1000
particles to pass the first criterion outlined above. In this case, using a
minimal group length of order 500 is beneficial over the more traditional
value of 32 as it will reduce the number of haloes to consider by about a
factor 10 (assuming a normal cosmology) whilst still being far enough from
the actual user-defined mass limit to be safe.

On-the-fly Friends-Of-Friends for snapshot output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It may be desirable to include FOF group membership information for each
particle in the output snapshots even when black hole seeding is not in use.
This can be achieved by setting the ``invoke_fof`` parameter in the 
``Snapshots`` section of the parameter file.

FOF will be run just before each snapshot is written and the snapshot will
include a dataset which specifies which FOF group each particle belongs to.
