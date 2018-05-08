Canonical Moon-Forming Giant Impact
===================================

A version of the canonical moon-forming giant impact of Theia onto the early
Earth (Canup 2004; Barr 2016). Both bodies made of Tillotson iron and granite.


Code Setup
----------

In `swiftsim/`:

`$ ./configure --with-hydro=minimal-multi-mat --with-equation-of-state=tillotson \ `
`   --disable-vec --with-kernel=wendland-C2`
`$ make`

In `swiftsim/examples/MoonFormingImpact/`:

`$ ./run.sh`


To do
-----

* Get it working!
* Add analysis code
* Make initial conditions available
