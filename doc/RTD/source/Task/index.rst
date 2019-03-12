.. Task
   Loic Hausammann 17th July 2018

.. _task:
   
Task System
===========

This section of the documentation includes information on the task system
available in SWIFT, as well as how to implement your own task.

SWIFT can produce a graph containing all the dependencies using graphviz.
At the beginning of each simulation a ``csv`` file is generated and can be transformed into a ``png`` with the script ``tools/plot_task_dependencies.py``.
This script has also the possibility to generate a list of function calls for each task with the option ``--with-calls``.
You can convert the ``dot`` file into a ``png`` with the following command
``dot -Tpng dependency_graph.dot -o dependency_graph.png`` or directly read it with the python module ``xdot`` with ``python -m xdot dependency_graph.dot``.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   adding_your_own
