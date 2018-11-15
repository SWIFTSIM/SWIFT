.. Engine
   November 2018
   Mladen Ivkovic



.. _engine:

Engine
------------------

Can be found in ``swiftsim/src/engine.c``


Main purpose: Handle everything about advancing the simulation in time.
Think of it as a "SWIFT blob" that contains and runs the heart of SWIFT.




Calling sequence in ``swiftsim/examples/main.c``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


- Define ``struct engine e`` that contains all the stuff the engine needs to advance the simulation in time, like the number of threads on which to run on, the space with wich the runner is associated, it's threads and policy, a task scheduler, starting and ending times, minimal and maximal timesteps, the MPI rank...
- :term:`pin the main thread <pinning a thread>`: Pins an engine (for each MPI rank) to a processor using pthreads black magic.
- If restarting: The previously used engine struct will be read in again, then configured. (The configuration in detail is a nightmare, so unless you're explicitly tinkering with the engine, don't open this Pandora's box.)
- Else: 

    - Construct :ref:`engine policy <engine_policy>`, initialize it and then construct the engine itself. 
    - :ref:`Split the domain <domain_decomposition>`:

        - engine_split:
            - Do initial partition
            - create :term:`proxies <proxy>`
            - re-allocate and re-:term:`link <links>` particles to add buffer space in the arrays

        - engine_redistribute:
            - Figure out how many particles of which :term:`type <particle types>` each node needs to send where
            - Communicate this information to other nodes.
            - Exchange particles (asynchronously)
            - re-link the particle types

    - initialize particles:
        - set particles to a valid state by calling ``engine_first_init_particles``

            - test

        - construct all cells and tasks, get ready to start running by calling ``engine_rebuild``

            - todo

        - compute densities? ``engine_skip_force_and_kick``
        - initialize particle data by calling ``space_init_(g/s)parts``
        - launch engine
        - prepare all tasks for a new round, then call ``engine_skip_drift``
        - initialize particle data by calling ``space_init_(g/s)parts``
        - recover the integer end of the next time step: ``engine_collect_end_of_step``
        - do some checks on particles









.. _engine_policy:

Engine Policy
~~~~~~~~~~~~~~~~~~~~~~~~

text
