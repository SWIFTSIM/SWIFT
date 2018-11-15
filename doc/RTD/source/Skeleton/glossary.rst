.. Glossary for code skeleton
   November 2018
   Mladen Ivkovic


.. _glossary:

Glossary
-----------------

.. glossary::

    Affinity
        By default, when running multiple tasks on multicore architectures, tasks are allowed to run on any core, and also switch on which core they are run on.
        To boost performance, instead you can pin a task on a specific core, where the data the task will require is stored most closely in memory, thus avoiding unnecessary data transfers from other cores.

    inhibited particle
        TODO

    links
        Each particle is split into "2 pieces": ``part`` + ``gpart``. (For stars, it's ``gpart`` + ``spart``. For dark matter only, no hydro/stars are needed, so it's a standalone ``gpart``). This is done so you can split hydro/gravity calculations without working on the same part of memory while computing them simultaneously. However, the particle parts need to be linked properly, which is what is meant by a link.

    output parameter file
        contains a list of what particle properties exist that can be written down in snapshot. 
        See :ref:`Output_list_label`

    pinning a thread
        Assigning a thread to a specific core to run on. 
        See :term:`Affinity`

    proxy
        structs that represent other MPI ranks; They contain the data on what information needs to be sent to other ranks

    particle types
        multiple particle types are employed in SWIFT. This allows tasks to compute multiple kinds of interactions on particle simultaneously. There are  "normal" SPH particles, star particles and gravity particles. For more details, refer to the :ref:`particles` section.

    runner
        the threads. Structure that has the main properties of threads like ID, type, etc. 
        
    stop file
        Undocumented feature at the moment. A way to stop the simulation in a clean way: SWIFT will periodically check whether such a stop file exists, and stop the code cleanly when it does such that it can be restarted cleanly.
