Dam Break (Tillotson water, closed tank, 3D)
--------------

A hydrostatic column of water, 0.5 m wide and 1 m high, rests against the
`x_min` wall of a closed 2 m x 2 m tank and collapses under constant gravity.
The tank is a thin slab along y (0.25 m), making the flow quasi-2D. All six faces of
the box are reflective walls (`BoundaryConditions:type: reflective`), so no
particle should ever leave the box.

The water uses the Tillotson equation of state (material ID 102), which has
zero pressure at its reference density and hence supports a free surface.
Gravity is scaled up to 5000 m/s^2 so that the flow stays weakly compressible
(Mach ~0.1) without the tiny time-steps the sound speed of water would
otherwise impose. Times scale with the free-fall time sqrt(H / g) ~ 0.014 s.

The walls only mirror the particles that cross them. There are no ghost
particles, so the kernels of the particles close to a wall are truncated.

The number of particles is set when generating the initial conditions and
defaults to one million: `python3 makeIC.py --n_part 1000000`.

`python3 makeMovie.py` renders every snapshot in the directory into
`dam_break.mp4` (requires swiftsimio and ffmpeg). Lower `Snapshots:delta_time`
for a smoother movie.

Configuration options:

REMIX:
`--with-hydro=remix --with-equation-of-state=planetary --with-kernel=wendland-C2 --with-ext-potential=constant`
