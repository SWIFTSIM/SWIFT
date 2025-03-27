Kelvin--Helmholtz Instabilty (Earth-like, equal mass, 3D)
--------------

This is a 3D version of the Kelvin--Helmholtz instability. These initial
conditions are those used to produce the simulations presented by
Sandnes et al. (2025), section 4.4.

This test uses particles of equal mass and has sharp density and velocity
discontinuities. Equations of state and conditions are representative of those
within Earth's interior.

Note that the default resolution_eta parameter is consistent with the use of a
Wendland C2 kernel with ~100 neighbours.

Some examples of configuration options with different hydro schemes:

REMIX:
`--with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2`

"Traditional" SPH (tSPH):
`--with-hydro=planetary --with-equation-of-state=planetary --with-kernel=wendland-C2`
