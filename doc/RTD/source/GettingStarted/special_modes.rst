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

Gravity glasses
~~~~~~~~~~~~~~~

For many problems in cosmology, it is important to start a simulation with no
initial noise in the particle distribution. Such a "glass" can be created by
starting from a random distribution of particles and running with the sign of
gravity reversed until all the particles reach a steady state. To run SWIFT in
this mode, configure the code with ``--enable-glass-making``.

Note that this will *not* generate the initial random distribution of
particles. An initial condition file with particles still has to be provided.

