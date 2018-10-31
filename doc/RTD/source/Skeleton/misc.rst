.. Misc stuff that I don't know where else to pu
   Mladen Ivkovic


Miscellaneous Stuff
-------------------------

Here's a list of stuff that I didn't know where else to put, but seem to be important concepts.
In alphabetical order.



.. glossary::

    Affinity
        By default, when running multiple tasks on multicore architectures, tasks are allowed to run on any core, and also switch on which core they are run on.
        To boost performance, instead you can pin a task on a specific core, where the data the task will require is stored most closely in memory, thus avoiding unnecessary data transfers from other cores.


    Pinning a thread
        Assigning a thread to a specific core to run on. See :term:`Affinity`
