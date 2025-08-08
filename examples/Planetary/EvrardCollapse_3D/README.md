Evrard Collapse 3D (Planetary)
===============

This is a copy of `/examples/HydroTests/EvrardCollapse_3D` for testing the 
Planetary and REMIX hydro schemes with the planetary ideal gas equation of state.

The Planetary hydro scheme results should be highly similar to the Minimal
hydro scheme, though  slightly different because of the higher viscosity beta
used here. To recover  the Minimal scheme behaviour, edit `const_viscosity_beta`
from 4 to 3 in `src/hydro/Planetary/hydro_parameters.h`.

Note that the default resolution_eta parameter is consistent with the use of a
Wendland C2 kernel with ~100 neighbours.

Some examples of configuration options with different hydro schemes:

REMIX:
`--with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2`

"Traditional" SPH (tSPH):
`--with-hydro=planetary --with-equation-of-state=planetary --with-kernel=wendland-C2`
