.. Misc stuff that I don't know where else to pu
   Mladen Ivkovic


.. _misc_glossary:

Miscellaneous Stuff and Glossary
--------------------------------------

Here's a list of stuff that I didn't know where else to put, but seem to be important concepts.
In alphabetical order.



.. glossary::

    Affinity
        By default, when running multiple tasks on multicore architectures, tasks are allowed to run on any core, and also switch on which core they are run on.
        To boost performance, instead you can pin a task on a specific core, where the data the task will require is stored most closely in memory, thus avoiding unnecessary data transfers from other cores.

    Links
        Each particle is split into "2 pieces": ``part`` + ``gpart``. (For stars, it's ``gpart`` + ``spart``. For dark matter only, no hydro/stars are needed, so it's a standalone ``gpart``). This is done so you can split hydro/gravity calculations without working on the same part of memory while computing them simultaneously. However, the particle parts need to be linked properly, which is what is meant by a link.

    output parameter file
        contains a list of what particle properties exist that can be written down in snapshot. See :ref:`Output_list_label`

    pinning a thread
        Assigning a thread to a specific core to run on. See :term:`Affinity`

    runner
        the threads. Structure that has the main properties of threads like ID, type, etc. 
        

    stop file
        Undocumented feature at the moment. A way to stop the simulation in a clean way: SWIFT will periodically check whether such a stop file exists, and stop the code cleanly when it does such that it can be restarted cleanly.
