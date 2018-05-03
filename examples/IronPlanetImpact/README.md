Iron Earths Giant Impact
========================

This is a first test of the planetary equations of state, colliding two Earth-
mass planets made of Tillotson iron.


Code Setup
----------

In `swiftsim/`:

`$ ./configure --with-hydro=minimal-multi-mat --with-equation-of-state=tillotson --disable-vec`
`$ make`

In `swiftsim/examples/IronPlanetImpact/`:

`$ ./run.sh`


To do
-----

* Get it working!
* Add analysis code
* Make initial conditions available
