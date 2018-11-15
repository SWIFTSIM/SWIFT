.. particle initialisation
   November 2018
   Mladen Ivkovic




.. _particles:

Particles
-----------------------------


General Information
~~~~~~~~~~~~~~~~~~~~~~~~~


**Particle Types**


SWIFT utilizes multiple particle types, defined as different structs. The exact definition used depends on the methods you chose to use. The following are descriptions of the defaut choice: 'minimal' hydro and 'default' gravity.
For details on others, have a look at the files themselves. The particle types are and they contain: (NOTE: 'ALFS' is an abbreviation for 'At Last Full Step')

- ``xpart``: a particle type containing particle fields not needed during the SPH loops over neighbours.

  .. code :: text

    offset between current position and position at last tree requild
    offset between the current position and position at the last sort
    velocity ALFS
    gravitational acceleration ALFS
    internal energy ALFS
    additional data for cooling information

- ``part``: hydro (SPH) particles
    
  .. code-block :: none

    particle ID
    pointer to corresponding gravity particle
    position
    velocity
    acceleration
    mass
    smoothing length
    internal energy
    time derivative of internal energy
    density
    force and density related quantities in a struct, stored as an union
    chemistry information in a chemistry_part_data struct [if no chemistry is selected, this will use no memory]
    time bin / time step length of particle

- ``gpart``: gravity particles

  .. code-block :: none
    
    particle ID or array index of hydro particle with which this gpart is linked
    position
    velocity
    acceleration
    mass
    time bin/time step length of particle
    type of gravity particle (dark matter, gas, star...)

- ``spart``: star particles

  .. code-block :: none
    
    particle ID
    pointer to corresponding gravity particle
    particle position
    offset between current position and position at last tree rebuild
    particle velocity
    mass
    cutoff radius
    time bin/time step length of particle
    struct with density related quantities


The multiple particle types are used to separate various interactions and enables the tasks to compute the interactions independently of each other simultaneously on the same particle.
Hydro and star particles have a corresponding gravity particle, respectively, which are linked by pointers.
As you can see, some particle properties are present multiple times in different particle types. While requiring more memory than absolutely necessary, this allows to boost performance by keeping memory local to each CPU as well as easier MPI communications because particles of each type encapsulate the data they need in their own structs.










.. _particle_init:

Particle Initialisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section describes how particles are stored and handled during initialisation, and what state they're left in before the main simulation loop starts, in case the simulation is started from scratch or restarted.


General Information
~~~~~~~~~~~~~~~~~~~~~


- the pointers to the ``gparts``, ``parts`` and ``sparts`` arrays are defined at the beginning of ``main.c`` to make sure that they're always in scope. They are pointers to arrays of the respective types and are initialized to ``NULL``.
- 




Fresh Start
~~~~~~~~~~~~~~~~~~~~
- the particles are read in in ``main.c`` after all necessary stuff has been set up. Depending on what libraries you use, they will be read in slightly differently:

    - If you have parallel HDF5: ``read_ic_parallel`` from ``parallel_io.c`` is called

        - Every MPI task reads in the total number of particles of each type
        - Particles are assigned to a MPI task by evenly dividing up contiguous sections of every particle type among all MPI tasks. Every MPI task now knows that it will initially contain ``N[particle_type]`` particles of type ``particle type``, which will be a fraction of ``N_total[particle_type]``.
        - particle arrays for ``parts``, ``gparts`` and ``sparts`` are allocated with size ``N[particle_type]`` and set to zero.
        - prepare which particle fields to read from the IC files by calling ``hydro_read_particles``, ``darkmatter_read_particles`` and ``stars_read_particles``
        - every MPI task reads in it's own share of particle data using HFD5 magic.

        TODO: left off at parallel_io.c line 881. Check how exactly arrays are populated.

    - If you don't have parallel HDF5, but use MPI:

    - If you only use non-parallel HDF5:





Restart
~~~~~~~~~~~~~~~~~~~~
