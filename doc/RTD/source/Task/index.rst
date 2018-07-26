.. Task
   Loic Hausammann 17th July 2018

.. _task:
   
Task System
===========

This section of the documentation includes information on the task system
available in SWIFT, as well as how to implement your own task.

SWIFT produces at the beginning of each simulation a ``dot`` file (see the graphviz library for more information).
It contains the full hierarchy of tasks used in this simulation.
You can convert the ``dot`` file into a ``png`` with the following command
``dot -Tpng dependency_graph.dot -o dependency_graph.png``.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   adding_your_own
