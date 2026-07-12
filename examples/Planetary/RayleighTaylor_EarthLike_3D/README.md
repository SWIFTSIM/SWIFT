Rayleigh--Taylor Instabilty (Earth-like, equal mass, 3D)
--------------

This is a 3D version of the Rayleigh--Taylor instability. These initial
conditions are those used to produce the simulations presented by
Sandnes et al. (2025), section 4.6.

This test uses particles of equal mass and has a sharp density. Equations
of state and conditions are representative of those within Earth's interior.

Note that the default resolution_eta parameter is consistent with the use of a
Wendland C2 kernel with ~100 neighbours.

Some examples of configuration options with different hydro schemes:

REMIX:
`--with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2 --with-ext-potential=constant --with-forcing=boundary-particles`

"Traditional" SPH (tSPH):
`--with-hydro=planetary --with-equation-of-state=planetary --with-kernel=wendland-C2 --with-ext-potential=constant --with-forcing=boundary-particles`

The flag `--with-forcing=boundary-particles` requires that the maximum boundary particle ID is
specified in the parameter file, for example:

BoundaryParticles:
  boundary_particle_max_id:        20992

This value is provided by the script makeIC.py and depends on the resolution.
