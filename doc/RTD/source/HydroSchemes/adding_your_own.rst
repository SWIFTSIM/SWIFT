.. Adding Hydro Schemes
   Josh Borrow, 5th April 2018


Adding Hydro Schemes
====================

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents:

SWIFT is engineered to enable you to add your own hydrodynamics schemes easily.
We enable this through the use of header files to encapsulate each scheme.

Note that it's unlikely you will ever have to consider parallelism or 'loops over
neighbours' for SWIFT; all of this is handled by the tasking system. All we ask
for is the interaction functions that tell us how to a) compute the density
and b) compute forces.


Getting Started
---------------

The hydro schemes are stored in ``src/hydro``. You will need to create a folder
with a sensible name that you are going to store your scheme-specific information
in. Then, you will need to create the following files:

+ ``hydro.h``, which includes functions that are applied to particles at the end
  of the density loop and beginning of the force loop, along with helper functions
+ ``hydro_debug.h``, which includes a quick function that prints out your particle
  properties for debugging purposes
+ ``hydro_iact.h`` that includes the interaction functions
+ ``hydro_io.h`` which includes the information on what should be read from the
  initial conditions file, as well as written to the output files
+ ``hydro_part.h`` which includes your particle definition. SWIFT uses an array-of
  -structures scheme.


``hydro.h``
-----------

As previously noted, ``hydro.h`` includes the helper functions for your scheme. You
will need to 'fill out' the following:

+ ``hydro_get_comoving_internal_energy(p)`` which returns the comoving internal energy
  of your particles (typically this will just be ``p->u``).
+ ``hydro_get_physical_internal_energy(p, cosmo)`` which returns the physical internal
  energy. You can use the ``a_factor_internal_energy`` from the ``cosmology`` struct.
+ ``hydro_get_comoving_pressure(p)`` which returns the comoving pressure.
+ ``hydro_get_comoving_entropy(p)`` which returns the comoving entropy.
+ ``hydro_get_physical_entropy(p, cosmo)`` which returns the physical entropy. In our
  formalism, usually there is no conversion factor here so it is the same as the
  comoving version.
+ ``hydro_get_comoving_soundspeed(p)`` which returns the comoving sound speed.
+ ``hydro_get_physical_soundspeed(p, cosmo)`` which returns the physical sound
  speed. You can use the ``a_factor_sound_speed``.
+ ``hydro_get_comoving_density(p)`` which returns the comoving density.
+ ``hydro_get_physical_density(p, cosmo)`` which returns the physical density.
  You can use the ``a3_inv`` member of the ``cosmology`` struct.
+ ``hydro_get_mass(p)`` returns the mass of particle ``p``.
+ ``hydro_get_drifted_velocities(p, xp, dt_kick_hydro, dt_kick_grav, v[3])`` gets
  the drifted velocities; this is just ``a_hydro * dt_kick_hydro`` + ``a_grav *
  dt_kick_grav`` in most implementations.
+ ``hydro_get_energy_dt(p)`` returns the time derivative of the (comoving) internal
  energy of the particle.
+ ``hydro_set_energy_dt(p)`` sets the time derivative of the (comoving) internal
  energy of the particle.
+ ``hydro_compute_timestep(p, xp, hydro_props, cosmo)`` returns the timestep for 
  the hydrodynamics particles.
+ ``hydro_timestep_extra(p, dt)`` does some extra hydro operations once the
  physical timestep for the particle is known.
+ ``hydro_init_part(p, hydro_space)`` initialises the particle in preparation for
  the density calculation. This essentially sets properties, such as the density,
  to zero.
+ ``hydro_end_density(p, cosmo)`` performs operations directly after the density
  loop on each particle. Note that you will have to add a particle's self-contribution
  at this stage as particles are never 'interacted' with themselves.
+ ``hydro_part_has_no_neighbours(p, xp, cosmo)`` resets properties to a sensible
  value if a particle is found to have no neighbours.
+ ``hydro_prepare_force(p, xp, cosmo)`` is computed for each particle before the
  force loop. You can use this to pre-compute particle properties that are used
  in the force loop, but only depend on the particle itself.
+ ``hydro_reset_acceleration(p)`` resets the acceleration variables of the particles
  to zero in preparation for the force loop.
+ ``hydro_predict_extra(p, xp, dt_drift, dt_therm)`` predicts extra particle properties
  when drifting, such as the smoothing length.
+ ``hydro_end_force(p, cosmo)`` is called after the force loop for each particle and 
  can be used to e.g. include overall factors of the smoothing length.
+ ``hydro_kick_extra(p, xp, dt_therm)`` kicks extra variables.
+ ``hydro_convert_quantities(p, xp)`` converts quantities at the start of a run (e.g.
  internal energy to entropy).
+ ``hydro_first_init_part(p, xp)`` is called for every particle at the start of a run
  and is used to initialise variables.


``hydro_debug.h``
-----------------

TBD


``hydro_iact.h``
----------------

TBD


``hydro_io.h``
--------------

TBD


``hydro_part.h``
----------------

TBD

