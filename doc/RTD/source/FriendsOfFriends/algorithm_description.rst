.. Friends Of Friends
   Matthieu Schaller 15th June 2019

.. _fof_algorithm_description_label:

Friends-Of-Friends Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Friends-Of-Friends (FOF) is a common tool used to identify groups of
particles in a simulation. This can, for instance, be used to identify
haloes in a cosmological simulation.

In practice, this is done by using a *linking length* ``l`` and
demanding that any particle that finds another particle within a
distance ``l`` is linked to it to form a group. A particle is linked
directly to all other particles within a distance ``l`` (its
*friends*) and indirectly to all particles that are linked to its
friends (its *friends-of-friends*). This creates networks of linked particles
which are called *groups*. The size (or length) of
a group is the number of particles in that group. If a particle does not
find any other particle within ``l`` then it forms its own group of
size 1. **For a given distribution of particles the resulting list of
groups is unique and unambiguously defined.**

In our implementation, we use three separate categories influencing their
behaviour in the FOF code:

- ``linkable`` particles which behave as described above.
- ``attachable`` particles which can `only` form a link with the `nearest` ``linkable`` particle they find.
- And the others which are ignored entirely.

The category of each particle type is specified at run time in the parameter
file. The classic scenario for the two categories is to run FOF on the dark
matter particles (i.e. they are `linkable`) and then attach the gas, stars and
black holes to their nearest DM (i.e. the baryons are `attachable`).

Small groups are typically discarded, the final catalogue only contains
objects with a length above a minimal threshold, typically of the
order of ``20`` particles. Smaller groups can often be spurious.

Once the groups have been identified, properties can be computed for
each of them. The total mass or the centre of mass are common
examples. These are then stored in catalogues alongside a unique
identifier for each group.

SWIFT implements FOF using a Union-Find approach. It also exploits the
domain decomposition and tree structure that is created for the other
parts of the code. The tree can be easily used to find neighbours of
particles within the linking length.

Depending on the application, the choice of linking length and minimal group
size can vary. For cosmological applications, bound structures (dark matter
haloes) are traditionally identified using a linking length expressed as
:math:`0.2` of the mean inter-particle separation :math:`d` in the simulation
which is given by :math:`d = \sqrt[3]{\frac{V}{N}}`, where :math:`N` is the
number of particles in the simulation and :math:`V` is the simulation
(co-moving) volume. Experience shows that this produces groups that are similar
to the commonly adopted (but much more complex) definition of virialised
haloes. A minimal group length of :math:`32` is often adopted in order to get a
robust catalogue of haloes and compute a good halo mass function.  Usually only
dark matter particles are considered for the number :math:`N`. In practice, the
mean inter-particle separation is evaluated based on the cosmology adopted in
the simulation.  We use: :math:`d=\sqrt[3]{\frac{m_{\rm DM}}{\Omega_{\rm cdm}
\rho_{\rm crit}}}` for simulations with baryonic particles and
:math:`d=\sqrt[3]{\frac{m_{\rm DM}}{(\Omega_{\rm cdm} + \Omega_{\rm b})
\rho_{\rm crit}}}` for DMO simulations. In both cases, :math:`m_{\rm DM}` is the
mean mass of the DM particles. Using this definition (rather than basing in on
:math:`N`) makes the code robust to zoom-in scenarios where the entire volume is
not filled with particles.

For non-cosmological applications of the FOF algorithm, the choice of
the linking length is more difficult and left to the user. The choice
of the minimal group size to consider is also application dependent.

In SWIFT, FOF is also used to identify haloes in which to seed black
holes in a cosmological run. When used in this mode, the code does not
produce any outputs but uses the catalogue internally to convert some
gas particles into black hole particles.

