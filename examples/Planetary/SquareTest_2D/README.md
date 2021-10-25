Square Test 2D (Planetary)
==============

This is a copy of `/examples/HydroTests/SquareTest_2D` for testing the Planetary 
hydro scheme with the planetary ideal gas equation of state. 

This is a very challenging test that aims to figure out
if contact discontinuities are properly handled. If there
is residual surface tension, then the square will quickly
become a sphere. Otherwise, it will remain a square. For
more information see Hopkins' 2013 and 2015 papers.

There are two initial condition generation files present.
For the SWIFT method of finding an un-mass weighted number
of particles in the kernel, it makes more sense to have
different mass particles (makeICDifferentMasses.py). For
comparison to previous methods, we also provide a script
that creates initial conditions with a different density
of particles, all with equal masses, in the square and
outside of the square.

If you do not have the swiftsimio library, you can use
the plotSolutionLegacy.py to plot the solution.

The results should be highly similar to the Minimal hydro scheme, though 
slightly different because of the higher viscosity beta used here. To recover 
the Minimal scheme behaviour, edit `const_viscosity_beta` from 4 to 3 in 
`src/hydro/Planetary/hydro_parameters.h`.

This requires the code to be configured to use the planetary hydrodynamics 
scheme and equations of state: 
`--with-hydro=planetary --with-equation-of-state=planetary --with-hydro-dimension=2`
