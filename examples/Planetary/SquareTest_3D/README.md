Square Test 3D
--------------

This is a 3D version of the "square test" consisting of a
cube of high-density material in pressure equilibrium with
a surrounding low-density region. These initial conditions
are thoe used to produce the simulations presented by
Sandnes et al. (2025).

This test is used to investigate spurious surface tension-like 
effects from sharp discontinuities in a system that should be 
in static equilibrium. There are two initial condition generation
files to test both an equal spacing scenario, i.e., with 
different particle masses in the two regions, and an equal mass 
scenario. The significant contributions from both smoothing and 
discretisation error at the density discontinuity makes the equal 
mass test particularly challenging for SPH.

Note that the default resolution_eta parameter is consistent
with the use of a Wendland C2 kernel with ~100 neighbours.

Some examples of configuration options with different hydro schemes: 

REMIX:
'--with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2'

'Traditional SPH (tSPH)' in Sandnes et al. (2025): 
'--with-hydro=planetary --with-equation-of-state=planetary --with-kernel=wendland-C2'

SPHENIX: 
'./configure --with-hydro=sphenix --with-kernel=wendland-C2'
