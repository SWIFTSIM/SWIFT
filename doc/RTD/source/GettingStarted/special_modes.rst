.. Special modes
   Matthieu Schaller, 20/08/2018

Special modes
=============

SWIFT comes with a few special modes of operating to perform additional tasks.

Static particles
~~~~~~~~~~~~~~~~

For some test problems it is convenient to have a set of particles that do not
perceive any gravitational forces and just act as sources for the force
calculation. This can be achieved by configuring SWIFT with the option
``--enable-no-gravity-below-id=N``. This will zero the *accelerations* of all
particles with ``id`` (strictly) lower than ``N`` at every time-step. Note that
if these particles have an initial velocity they will keep moving at that
speed.

This will also naturally set their time-step to the maximal value
(``TimeIntegration:dt_max``) set in the parameter file.

A typical use-case for this feature is to study the evolution of one particle
orbiting a static halo made of particles. This can be used to assess the
quality of the gravity tree and time integration. As more particles are added
to the halo, the orbits will get closer to the analytic solution as the noise
in the sampling of the halo is reduced.

Note also that this does not affect the hydrodynamic forces. This mode is
purely designed for gravity-only accuracy tests.

Besides the pure gravity mode, swift also has the boundary particle mode,
this mode turns off both the gravity forces and hydro forces for all
particles. Because gas particles only receive hydro this mode only impacts
gas particles more strictly than other particles. This mode can be
activated using ``--enable-boundary-particles=N``. This will zero the
gravitational and hydrodynamic *accelerations* of all particles with ``id``
(strictly) lower than ``N`` at every time-step. Still if particles have an
initial velocity they will keep moving at that speed. This compilation
option also activates ``--enable-no-gravity-below-id=N``. 

Typical use cases are hydrodynamical experiments that have boundaries. 

Both options discussed above only cancel accelerations and have no impact
on what can happen in any code module that directly changes the velocity of
the boundary or no gravity particles. An example of this is momentum
injection of stellar winds for example. If we additionally want to keep the
boundary particles fixed at the same position for the whole simulation we can
use the ``--enable-fixed-boundary-particles=N`` compile option, this option
explicitly sets the velocity of the boundary particles to zero every time
step on top of also zeroing the accelerations.

A typical use cases is an idealized setup in which we have a black hole in
the centre of a galaxy that accretes gas but is not allowed to move from
the momentum recoil of the gas it swallows.

Gravity glasses
~~~~~~~~~~~~~~~

For many problems in cosmology, it is important to start a simulation with no
initial noise in the particle distribution. Such a "glass" can be created by
starting from a random distribution of particles and running with the sign of
gravity reversed until all the particles reach a steady state. To run SWIFT in
this mode, configure the code with ``--enable-glass-making``.

Note that this will *not* generate the initial random distribution of
particles. An initial condition file with particles still has to be provided.

Gravity force checks
~~~~~~~~~~~~~~~~~~~~

To test the accuracy of the gravitational forces approximated by the code,
SWIFT can be configured with the option to additionally compute the "exact"
forces for each active particle during each timestep. Here the exact forces are
simply the Newtonian sum, i.e.,

:math:`\vec{F}_{i,j} = \sum^{n}_{i \neq j} \frac{G m_i m_j}{\vec{r}_{i,j}^2}.`

To run SWIFT in this mode, configure the code with
``--enable-gravity-force-checks=N``, which means that the exact forces will be
computed for every :math:`N^{th}` particle based on their ID (i.e., to compute
the exact forces for all particles set ``N=1``).

Two `.dat` files will be output during each timestep, one containing the forces
(really it is the accelerations that are stored) as computed by ``_swift_``, and
another containing the ``_exact_`` forces. The total force (``a_swift``), along
with the contributed force from each component of the tree (P2P, M2P and M2L)
and the FFT mesh if periodic (PM) is stored (i.e., ``a_swift[0]`` = ``a_p2p[0]`` +
``a_m2p[0]`` + ``a_m2l[0]`` + ``a_PM[0]``, for the :math:`x` component). In addition,
the number of particles contributing to each force component is also stored
(these numbers will add up to :math:`n-1`).   

This mode will slow down the code *considerably*, and it is not recommended to
be run on problems with more than :math:`10^{5}` particles when
``--enable-gravity-force-checks=1``. For larger runs, sampling a sub-set of
particles via the argument ``N`` of the configuration option is recommended.
This mode must be run on a single node/rank, and is primarily designed for pure
gravity tests (i.e., DMO).
