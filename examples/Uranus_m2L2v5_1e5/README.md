Uranus Giant Impact, low ang. mom.
==================================

As in Kegerreis et al. (2018): A relatively head-on impact onto the early Uranus, 
with an impactor mass of 2 Earth masses, angular momentum of 2e36 kg m^2 s^-1,
and velocity at infinity of 5 km s^-1.

Code Setup
----------

In `swiftsim/`:

`$ ./configure --with-hydro=minimal-multi-mat --with-equation-of-state=planetary \ `
`   --disable-vec --with-kernel=wendland-C2`
`$ make`

In `swiftsim/examples/Uranus_m2L2v5_1e5/`:

`$ ./run.sh`


To do
-----

* Get it working!
* Add analysis code
* Make initial conditions available
