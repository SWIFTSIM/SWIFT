.. AnalysisTools
   Loic Hausammann 20th March 2019

.. _analysistools:

Analysis Tools
==============

Task dependencies
-----------------

At the beginning of each simulation the file ``dependency_graph.csv`` is generated and can be transformed into a ``dot`` and a ``png`` file with the script ``tools/plot_task_dependencies.py``.
It requires the ``dot`` package that is available in the library graphviz.
This script has also the possibility to generate a list of function calls for each task with the option ``--with-calls`` (this list may be incomplete).
You can convert the ``dot`` file into a ``png`` with the following command
``dot -Tpng dependency_graph.dot -o dependency_graph.png`` or directly read it with the python module ``xdot`` with ``python -m xdot dependency_graph.dot``.


Cell graph
----------

An interactive graph of the cells is available with the configuration option ``--enable-cell-graph``.
During a run, SWIFT will generate a ``cell_hierarchy_*.csv`` file per MPI rank.
The command ``tools/make_cell_hierarchy.sh cell_hierarchy_*.csv`` merges the files together and generates the file ``cell_hierarchy.html``
that contains the graph and can be read with your favorite web browser.

With chrome, you cannot access the files directly, you will need to either access them through an existing server (e.g. public http provided by your university)
or install ``npm`` and then run the following commands

.. code-block:: bash
   
   npm install http-server -g
   http-server .

Now you can open the web page ``http://localhost:8080/cell_hierarchy.html``.
