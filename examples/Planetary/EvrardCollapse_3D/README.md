Evrard Collapse 3D (Planetary)
===============

This is a copy of `/examples/HydroTests/EvrardCollapse_3D` for testing the 
Planetary hydro scheme with the planetary ideal gas equation of state. 

The results should be highly similar to the Minimal hydro scheme, though 
slightly different because of the higher viscosity beta used here. To recover 
the Minimal scheme behaviour, edit `const_viscosity_beta` from 4 to 3 in 
`src/hydro/Planetary/hydro_parameters.h`.

This requires the code to be configured to use the planetary hydrodynamics 
scheme and equations of state: 
`--with-hydro=planetary --with-equation-of-state=planetary`
