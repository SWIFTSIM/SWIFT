.. Engine
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
- If restarting: The previously used engine struct will be read in again, then configured. (The configuration in detail is a nightmare, so unless you're explicitly tinkering with the engine, just skip it.)
- Else: 

    - Construct engine policy, initialize it and then construct the engine itself. 
    - (Later) split the domain

        - **engine_split**:
        - Do initial partition
        - create :term:`proxies <proxy>`
        - re-allocate and re-:term:`link <links>` particles to add buffer space in the arrays

        - **engine_redistribute**:
        - 



