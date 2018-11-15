.. main.c
   November 2018
   Mladen Ivkovic



Main Function
------------------

located at ``swiftsim/examples/main.c``

- Initialize MPI
- read command line arguments and options
- write :term:`output parameter file` if necessary (i.e. if it was given.)
- checking, announcing, and setting up stuff : 

    - is parameter file given
    - check whether you chose appropriate flags
    - set up some other minor stuff, :term:`pin the main thread <pinning a thread>`
    - make some announcements for the log/screen
    - read and broadcast parameter file
    - check whether you can write data
    - Initialize :ref:`domain decomposition <domain_decomposition>` (without actually doing it)
    - check where to read/write restart files
    - setup run parameters, i.e. how often to check for the :term:`stop file` , maximal wall-clock time, whether to resubmit the run...
- if restarting from some snapshot, check that restart files exist, read them and start the :ref:`engine`. Otherwise: 
    
    - look for initial conditions (IC), and read them
    - initialize units, constants, physics (cosmology, hydro, stars...)
    - initialize space, distribute particles across it
    - initialize further physics/chemistry stuff (ext. potential, cooling, chemistry, feedback...)
    - set up engine policies, initialize :ref:`engine`.

- dump parameters as used
- main simulation loop:

    - reset timers
    - take an engine step
    - check whether you need to stop
    - dump whatever needs to be dumped
    - clean up
    - rinse and repeat

