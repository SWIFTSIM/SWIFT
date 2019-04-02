.. AnalysisTools
   Loic Hausammann 20th March 2019
   Peter W. Draper 28th March 2019

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

Memory usage reports
--------------------

When SWIFT is configured using the ``--enable-memuse-reports`` flag it will
log any calls to allocate or free memory that make use of the
``swift_memalign()``, ``swift_malloc()``, ``swift_calloc()`` and
``swift_free()`` functions and will generate a report at the end of each
step. It will also attempt to dump the current memory use when SWIFT is
aborted by calling the ``error()`` function. Failed memory allocations will be
reported in these logs.

These functions should be used by developers when allocating significant
amounts of memory -- so don't use these for high frequency small allocations.
Each call to the ``swift_`` functions differs to the standard calls by the
inclusion of a "label", this should match between allocations and frees and
ideally should be a short label that describes the use of the memory, i.e.
"parts", "gparts", "hydro.sort" etc.

Calls to external libraries that make allocations you'd also like to log
can be made by calling the ``memuse_log_allocation()`` function directly.

The output files are called ``memuse_report-step<n>.dat`` or
``memuse_report-rank<m>-step<n>.dat`` if running using MPI. These have a line
for each allocation or free that records the time, memory address, step,
whether an allocation or free, the label and when an allocation, the amount of
memory. The comments in this file also record the actual memory use of the
process (including threads) as reported by the operating system at the end of
the step.

To post process these files into a memory used timeline and get a report of
the peak memory use, as well as the memory still in use at the end of the step
there is an basic analysis script ``analyse_memuse_logs.py`` and two wrappers
that process a directory of logs, these are ``./process_memuse_logs.sh`` and
``./process_memuse_logs_MPI.sh`` for non-MPI and MPI runs respectively.

Note that the process scripts process each step individually and also process
all the logs as a single sequence. When interpreting these some care should be
taken as they are not all the memory allocated, just important allocations in
SWIFT and when looking at a single step the context of any previous steps is
not used, so you only see allocations made in that step and the effect of any
matching frees (so allocations made in previous steps that are freed in this
step will not be understood and will be ignored, you need the global analysis
to understand that).

