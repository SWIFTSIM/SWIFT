.. Planetary Removed Particles
   Jacob Kegerreis, 27th November 2020

.. _planetary_removed_particles:
   
Removed Particles
=================

If a particle leaves a non-periodic simulation box then it is removed. 
Currently, the particle's information is simply printed to stdout in a csv 
format (see ``hydro_remove_part()`` in ``src/hydro/Planetary/hydro.h``).
This could be extracted for analysis e.g. using grep:
``grep '## Removed' -A 1 --no-group-separator output.txt > removed.txt``, 
then read e.g. with python:

.. code-block:: python

    import numpy as np
    id, x, y, z, vx, vy, vz, mass, u, P, rho, h, mat_id, time = np.loadtxt(
        "removed.txt", delimiter=",", unpack=True
    )
    pos = np.transpose([x, y, z])
    vel = np.transpose([vx, vy, vz])
    id = id.astype(int)
    mat_id = mat_id.astype(int)