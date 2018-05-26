Canonical Moon-Forming Giant Impact
===================================

NOTE: This doesn't really work because the EOS are different to Canup (2004) so
the impactor just glances then flies away!

A version of the canonical moon-forming giant impact of Theia onto the early
Earth (Canup 2004; Barr 2016). Both bodies are here made of a (Tillotson) iron
core and granite mantle. Only ~10,000 particles are used for a quick and crude
simulation.

Setup
-----

In `swiftsim/`:

`$ ./configure --with-hydro=minimal-multi-mat --with-equation-of-state=planetary`
`$ make`

In `swiftsim/examples/MoonFormingImpact/`:

`$ ./get_init_cond.sh`

Run
---

`$ ./run.sh`

Output
------

`$ python plot.py`
`$ mplayer anim.mpg`

