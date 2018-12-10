.. VELOCIraptor output
   Folkert Nobels 12th of October

VELOCIraptor Output
===================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents: 

In general VELOCIraptor outputs six files per snapshot, of which 2 files are
for unbound particles specifically.  In this part we will explain what is
inside the different files.

Catalog_groups file
-------------------

The first output file of VELOCIraptor is the ``.catalog_group`` file,
this file contains all the information that is group specific, and does not go
into depth of physical properties but only on numbers of particles and 
group sizes, the interesting data in the ``.catalog_group`` files are: 

+ The ``group_size``: gives a list of all the halos and the number of particles
  in the halo, this list is numbered from 0 until the number of groups minus
  one. It is important that the groups are not ordered in any way [#order]_.
  It is also important to note that the group size includes both the bound and
  unbound particles; always use the ``Offset`` and ``Offset_unbound`` data
  when reading from the ``catalog_particles`` files.
+ The ``Num_of_groups`` or ``Total_num_of_groups``: gives the total number of
  groups in the snapshot.
+ The ``Offset`` list: This list gives the offset off the particles. In the
  output of VELOCIraptor there is no file which has an ID for every particle
  and a corresponding group, rather the particles are ordered according to in
  which group they are. So if we want to access the particles in group 0, we
  need to look at the particles from ``Offset[0]`` until ``Offset[1]`` in the
  ``.catalog_particles`` hdf5 file. In general this means that for group N we
  need to look at particles ``Offset[N]`` until ``Offset[N+1]``. 
+ The ``Offset_unbound`` list: This list works exactly the same as the
  ``Offset`` list only this list is for the gravitational unbound particles.

Catalog_particles file
----------------------

The second file that is produced by VELOCIraptor is the ``.catalog_particles``
file, this file contains mainly all the IDs of the particles and has two
interesting parameters:

+ The ``Num_of_particles_in_groups`` and ``Total_num_of_particles_in_all_groups``
  parameter: Gives the total number of particles in the file or the total 
  number of particles that are in halos.
+ The ``Particle_IDs``: The list of particles as sorted by halo, in which halo
  the individual particles are present can be found by using the
  ``.catalog_group`` file and the corresponding ``Offset`` list. 

Besides the ``.catalog_particles`` file, there is also a
``.catalog_particles.unbound`` file, this file contains the same information
but only for the unbound particles, a particle can only be present in one of
these two lists. 

Extracting the particles in a given halo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``.catalog_particles`` file returns particle IDs that need to be matched
with those in your snapshot to find the particles in the file that you
wish to extract. The python snippet below should give you an idea of how to
go about doing this for the bound particles.

First, we need to extract the offset from the ``.catalog_group`` file, and
work out how many _bound_ particles are in our halo. We can do this by
looking at the next offset. Then, we can ID match those with the snapshot
file, and get the mask for the _positions_ in the file that correspond
to our bound particles. (Note this requires ``numpy > 1.15.0``).

.. code-block:: python
   :linenos:

   import numpy as np
   import h5py

   snapshot_file = h5py.File("swift_snapshot.hdf5", "r")
   group_file = h5py.File("velociraptor_output.catalog_group", "r")
   particles_file = h5py.File("velociraptor_output.catalog_particles", "r")

   halo = 100
   # Grab the start position in the particles file to read from
   halo_start_position = group_file["Offset"][halo]
   halo_end_position = group_file["Offset"][halo + 1]
   # We're done with that file now, best to close earlier rather than later
   group_file.close()

   # Get the relevant particle IDs for that halo; this includes particles
   # of _all_ types.
   particle_ids_in_halo = particles_file["Particle_IDs"][
       halo_start_position:halo_end_position
   ]
   # Again, we're done with that file.
   particles_file.close()

   # Now, the tricky bit. We need to create the correspondence between the
   # positions in the snapshot file, and the ids.

   # Let's look for the dark matter particles in that halo.
   particle_ids_from_snapshot = snapshot_file["PartType1/ParticleIDs"][...]

   _, indices_v, indices_p = np.intersect1d(
       particle_ids_in_halo,
       particle_ids_from_snapshot,
       assume_unique=True,
       return_indices=True,
   )

   # indices_p gives the positions in the particle file where we will find
   # the co-ordinates that we're looking for! To get the positions of all of
   # those particles,
   particle_positions_in_halo = snapshot_file["PartType1/Coordinates"][indices_p]


Catalog_parttypes file
----------------------

The third file that is produced by VELOCIraptor is the ``.catalog_parttypes``
file, this file contains the information what type of particle every particle
is, it is ordered the same as the ``Particle_IDs`` in ``.catalog_particles``. 
There are only two interesting parameters of the file which are:

+ The ``Num_of_particles_in_groups`` parameter: Gives the total number of
  particles in the file which are in a halo.
+ The ``Particle_types`` list: Gives a list of particles types similar to the
  snap shots (0 - gas, 1 - dm, 4 - stars).

Besides the ``.catalog_parttypes`` file, there is also a
``.catalog_parttypes.unbound`` file, this file contains this information for
the unbound particles.

Properties file
---------------

The fourth file is the ``.properties`` file, this file contains many physical
useful information of the corresponding halos. This can be divided in several
useful groups of physical parameters, on this page we have divided the several
variables which are present in the ``.properties`` file. This file has most 
physical interesting parameters of the halos.

Mass-Radius determination:
^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``.properties`` file contains many ways to determine the size and mass 
of the halos, in this subsection we will list several available variables in
the output of VELOCIraptor and we list several mass and radius parameters in
the output which are not classified as a mass-radius pair.

Critical Density related:
"""""""""""""""""""""""""

+ ``Mass_200crit``: The mass of a halo with an over density on average of
  :math:`\Delta=200` based on the critical density of the Universe 
  (:math:`M_{200}`).
+ ``R_200crit``: The :math:`R_{200}` radius of the halo based on the 
  critical density of the Universe

Mean Density related:
"""""""""""""""""""""

+ ``Mass_200mean``: The mass of a halo with an over density on average of
  :math:`\Delta=200` based on the mean density of the Universe 
  (:math:`M_{200}`).
+ ``R_200mean``: The :math:`R_{200}` radius of the halo based on the 
  mean density of the Universe.

Virial properties:
""""""""""""""""""

+ ``Mvir``: The virial mass of the halos.
+ ``Rvir``: The virial radius of the halo (:math:`R_{vir}`).

Bryan and Norman 1998 properties:
"""""""""""""""""""""""""""""""""

+ ``Mass_BN98``, The Bryan and Norman (1998) determination of the mass of the
  halo [#BN98]_. 
+ ``R_BN98``, the Bryan and Norman (1998) corresponding radius [#BN98]_.

Several Mass types:
"""""""""""""""""""
This is a list of masses which cannot be categorized as easy as the other 
properties.

+ ``Mass_FOF``: The friends-of-friends mass of the halos.
+ ``M_gas``: The gas mass in the halo.
+ ``Mass_tot``: The total mass of the halo
+ ``M_gas_30kpc``: The gas mass within 30 kpc of the halo centre.
+ ``M_gas_500c``: The gas mass of the over-density of 500 times the critical
  density
+ ``M_gas_Rvmax``: The gas mass within the maximum rotation velocity.

Several Radius types:
"""""""""""""""""""""

+ ``R_HalfMass``: Radius of half the mass of the halo.
+ ``R_HalfMass_gas``: Radius of half the gas mass of the halo.
+ ``R_size``:
+ ``Rmax``: 

Mass Structure of the Halos:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this subsection we listed the properties of the halos that are determining 
the mass structure of the halo, so the exact profile and the inertia tensor.

NFW profile properties:
"""""""""""""""""""""""
+ ``Xc``, ``Yc`` and ``Zc``: The x,y and z centre positions of the
  halos.
  
  Centres are calculated using first all particles belonging to the
  structure and then VELOCIraptor uses shrinking spheres to iterate to
  a centre, stopping once the sphere contains <10% of all the
  particles (this value can be changed to smaller amounts and there is
  also a minimum particle number which can also be changed).
  
+ ``Xc_gas``, ``Yc_gas``, ``Zc_gas``: The offset of the centre
  positions of the halo based on the gas, to find the position of the
  gas the offsets need to be added to ``Xc``, ``Yc`` and ``Zc``.

+ ``cNFW``: The concentration of the halo.

  This is calculated using Vmax and Vvir, not using a fitted profile.
  
+ ``VXc``, ``VYc`` and ``VZc`` are the velocities in the centre of the halo
  [#check]_.
+ ``VXc_gas``, ``VYc_gas`` and ``VZc_gas`` are the velocities of the gas  in
  the centre of the halo [#check]_.

Inertia Tensor properties:
"""""""""""""""""""""""""""

+ ``eig_ij``: Are the normalized eigenvectors of the inertia tensor.
+ The eigenvalue ratios: 

  1. ``q`` is the semi-major over major; 
  2. ``s`` is the minor over major.

+ ``eig_ij_gas``: Are the normalized eigenvectors of the inertia tensor for
  only the gas particles.
+ The eigenvalue ratios for only the gas, similar to all particles:

  1. ``q_gas`` is the semi-major over major for only gas; 
  2. ``s_gas`` is the minor over major for only gas.

Dynamical Structure of the Halos:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this subsection we list several properties that determine the dynamical
structure of the halo, like the angular momentum and the velocity dispersion
tensor.

Angular momentum and spin parameters:
"""""""""""""""""""""""""""""""""""""

+ ``lambda_b`` is the bullock spin parameter, see the paper by Bullock et al.
  (2001) [#Bullock]_. 
+ ``Lx``, ``Ly`` and ``Lz`` are the angular momentum of the halos, the 
  calculation includes all the particle types.
+ ``Lx_gas``, ``Ly_gas`` and ``Lz_gas`` are the angular momentum for only 
  the gas particles in the snapshot.

Velocity Dispersion related:
""""""""""""""""""""""""""""

+ The complete velocity dispersion tensor (:math:`\sigma_{ij}`) which has 
  an array per component which gives the value for all the halos. In 
  general these components are called ``veldisp_ij`` in which i and j are 
  given by ``x``, ``y`` or ``z``. This means that there are nine 
  components stored in the ``.properties`` file. This omits the fact 
  that the dispersion tensor by nature is a symmetric tensor. All the 
  components are given by: 
  ``veldisp_xx``, ``veldisp_xy``, ``veldisp_xz``, ``veldisp_yx``, 
  ``veldisp_yy``, ``veldisp_yz``, ``veldisp_zx``, ``veldisp_zy``, 
  and ``veldisp_zz`` [#velodisp]_.
+ ``sigV``, the scalar velocity dispersion which corresponds with the 
  trace of the velocity dispersion tensor 
  (:math:`\sigma = \text{Tr}(\sigma_{ij})`).


Energy properties of the halos:
"""""""""""""""""""""""""""""""

+ ``Ekin``, the kinetic energy of the halo.
+ ``Epot``, the potential energy of the halo.
+ ``Krot``, the rotational energy of the halo.
+ ``Krot_gas``, the rotational energy of the gas in the halo.


Halo and subhalo abstract variables:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this subsection we list the ID convention for subhalos and halos and 
some other abstract quantities of the halo which are not physical but 
rather properties of the simulations.

Structure types:
""""""""""""""""

+ ``ID`` is the halo ID.
+ ``Structuretype`` is the parameter that indicates what kind of structure 
  the current halo is. Halos have a structure type of ``10`` and subhalos
  have a structure type of ``15``.
+ ``hostHaloID``, indicates the halo ID number of the host halo, in the case
  that the halo has no parent (e.g. is the largest halo), the hostHaloID will
  be ``-1``.
+ ``numSubStruct``, the number of substructures or subhalos in the halo.

Particle types:
"""""""""""""""

+ ``npart`` is the number of particles in the halo (all types of particles).
+ ``n_gas`` is the number of gas particles in the halo.

Not specified parameters:
^^^^^^^^^^^^^^^^^^^^^^^^^

In this section we list parameters which cannot specifically be classified 
in a group.


Most Bound Particle (MBP):
""""""""""""""""""""""""""

+ ``ID_mbp``, the ID of the most bound particle in the halo.
+ ``Xcmbp``, ``Ycmbp`` and ``Zcmbp`` are the positions of the most bound 
  halo particle [#check]_.
+ ``VXcmbp``, ``VYcmbp`` and ``VZcmbp`` are the velocities of the most bound
  halo particle [#check]_.

.. [#order] In most cases more massive groups appear earlier in the list, but 
   this is not guaranteed for larger simulations. The order of the groups is 
   more a matter of the way that VELOCIraptor searches instead of a physical 
   reason.
.. [#center] This is not the average positions of the halos particles, but
   the halo position found by the VELOCIraptor algorithm. This includes a 
   fit for all the parameters including the gas particles or other types of
   particles.
.. [#velodisp] In the velocity dispersion tensor ( :math:`\sigma_{ij}` )  
   the following relations are satisfied between components:

   + :math:`\sigma_{xy}=\sigma_{yx}`
   + :math:`\sigma_{xz}=\sigma_{zx}`
   + :math:`\sigma_{yz}=\sigma_{yz}`
.. [#Bullock] The Bullock spin parameter is given by 
   :math:`\lambda = \frac{J}{\sqrt{2}MVR}`, for more information see 
   https://arxiv.org/abs/astro-ph/0011001. 
.. [#BN98] The Bryan and Norman (1998) paper can be found here: 
   https://arxiv.org/abs/astro-ph/9710107
.. [#check] Needs to be checked.
