Uranus Giant Impact
===================

A simple version of the low angular momentum impact onto the early Uranus shown
in Kegerreis et al. (2018), Fig. 2; with only ~10,000 particles for a quick and
crude simulation.

The collision of a 2 Earth mass impactor onto a proto-Uranus that can explain
the spin of the present-day planet, with an angular momentum of 2e36 kg m^2 s^-1
and velocity at infinity of 5 km s^-1 for a relatively head-on impact.

Both bodies have a rocky core and icy mantle, with a hydrogen-helium atmosphere
on the target as well. Although with this low number of particles it cannot be
modelled in any detail.

Setup
-----

In `swiftsim/`:

`$ ./configure --with-hydro=minimal-multi-mat --with-equation-of-state=planetary`

`$ make`

In `swiftsim/examples/UranusImpact/`:

`$ ./get_init_cond.sh`

Run
---

`$ ./run.sh`

Analysis
--------

`$ python plot.py`

`$ mplayer anim.mpg`

