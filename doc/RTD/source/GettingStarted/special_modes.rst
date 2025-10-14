.. Special modes
   Matthieu Schaller, 20/08/2018

Special modes
=============

SWIFT comes with a few special modes of operating to perform additional tasks.

Disabling particle types
~~~~~~~~~~~~~~~~~~~~~~~~

To save some meory, SWIFT can be compiled with support for reduced number of
particles. This will make many structures in the code discard the variables
related to these particles and hence lead to a leaner code. This can be useful,
for instance, for very large DMO runs.

This is achieved by compiling with one or more of these:
 * ``--with-hydro=none``
 * ``--with-stars=none``
 * ``--with-sinks=none``
 * ``--with-black-holes=none``

The code will then naturally complain if particles of the cancelled types are
found in the simulation.

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

SWIFT also includes functionality for boundary particles, which can be activated
by configuring with the flag ``--with-forcing=boundary-particles``. By default,
this disables both gravitational and hydrodynamic forces for particles with IDs
less than or equal to a specified maximum ID. This value, along with additional
parameters to enable alternative modes (the default options for which are listed
below), is set in the simulation parameter file using:

.. code-block:: yaml

   BoundaryParticles:
     boundary_particle_max_id:        N
     enable_fixed_position:           0
     enable_hydro_acceleration:       0
     enable_grav_acceleration:        0

where ``N`` would be replaced by an integer particle ID. With these options, if
boundary particles have an initial velocity, they will continue to move with
that velocity for the duration of the simulation.

Additional modes can be enabled:

- ``enable_fixed_position: 1`` sets boundary particle velocities to zero, in addition to
  accelerations, at each time step, fully fixing them in place, even if they had an
  initial velocity.
- ``enable_hydro_acceleration: 1`` allows hydrodynamic forces to act on boundary particles.
- ``enable_grav_acceleration: 1`` allows gravitational forces to act on boundary particles.

Typical use cases for boundary particles include experiments with fixed
boundaries and hydrodynamic test scenarios, such as the Rayleigh--Taylor instability
examples provided in SWIFT.

Note that certain physics modules may alter boundary particle velocities if they
update momenta outside of contributions from hydrodynamic or gravitational
accelerations. One example is the change in black hole momentum due to
accretion. The configuration option ``--enable-fixed-black-holes=N`` can be used
to ensure that black hole particles remain completely stationary. When enabled,
all black hole particles with IDs less than or equal to ``N`` are fully fixed in
position.

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

Two ``.dat`` files will be output during each timestep, one containing the forces
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
