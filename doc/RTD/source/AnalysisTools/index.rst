.. AnalysisTools
   Loic Hausammann 20th March 2019
   Peter W. Draper 28th March 2019

.. _Analysis_Tools:

Analysis Tools
==============

Task dependencies
-----------------

At the beginning of each simulation the file ``dependency_graph_0.csv`` is generated and can be transformed into a ``dot`` and a ``png`` file with the script ``tools/plot_task_dependencies.py``.
It requires the ``dot`` package that is available in the library graphviz.
This script has also the possibility to generate a list of function calls for each task with the option ``--with-calls`` (this list may be incomplete) and to describe at which level each task are run ``--with-levels`` (a larger simulation will provide more accurate levels).
You can convert the ``dot`` file into a ``png`` with the following command
``dot -Tpng dependency_graph.dot -o dependency_graph.png`` or directly read it with the python module ``xdot`` with ``python -m xdot dependency_graph.dot``.
If you wish to have more dependency graphs, you can use the parameter ``Scheduler:dependency_graph_frequency``. It defines how many steps are done in between two graphs.
While the initial graph is showing all the tasks/dependencies, the next ones are only showing the active tasks/dependencies.



Task levels
-----------------

At the beginning of each simulation the file ``task_level_0.txt`` is generated. 
It contains the counts of all tasks at all levels (depths) in the tree.
The depths and counts of the tasks can be plotted with the script ``tools/plot_task_levels.py``.
It will display the individual tasks at the x-axis, the number of each task at a given level on the y-axis, and the level is shown as the colour of the plotted point.
Additionally, the script can write out in brackets next to each tasks's name on the x-axis on how many different levels the task exists using the ``--count`` flag.
Finally, in some cases the counts for different levels of a task may be very close to each other and overlap on the plot, making them barely visible.
This can be alleviated by using the ``--displace`` flag: 
It will displace the plot points w.r.t. the y-axis in an attempt to make them better visible, however the counts won't be exact in that case.
If you wish to have more task level plots, you can use the parameter ``Scheduler:task_level_output_frequency``. 
It defines how many steps are done in between two task level output dumps.




Cell graph
----------

An interactive graph of the cells is available with the configuration option ``--enable-cell-graph``.
During a run, SWIFT will generate a ``cell_hierarchy_*.csv`` file per MPI rank at the frequency given by the parameter ``--cell-dumps=n``.
The command ``tools/make_cell_hierarchy.sh cell_hierarchy_0000_*.csv`` merges the files at time step 0 together and generates the file ``cell_hierarchy.html``
that contains the graph and can be read with your favorite web browser.

With most web browsers, you cannot access the files directly.
If it is the case, the cells will never appear (but everything else should be fine).
To solve this problem, you will need to either access them through an existing server (e.g. public http provided by your university)
or install ``npm`` and then run the following commands

.. code-block:: bash

   npm install http-server -g
   http-server .

Now you can open the web page ``http://localhost:8080/cell_hierarchy.html``.
When running a large simulation, the data loading may take a while (a few seconds for EAGLE_6).
Your browser should not be hanging, but will seems to be idle.

If you wish to add some information to the graph, you can do it by modifying the files ``src/space.c`` and ``tools/data/cell_hierarchy.html``.
In the first one, you will need to modify the calls to ``fprintf`` in the functions ``space_write_cell_hierarchy`` and ``space_write_cell``.
Here the code is simply writing CSV files containing all the required information about the cells.
In the second one, you will need to find the function ``mouseover`` and add the field that you have created.
You can also increase the size of the bubble through the style parameter ``height``.

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
for each allocation or free that records the time, step, whether an allocation
or free, the label, the amount of memory allocated or freed and the total of
all (labelled) memory in use at that time.

Comments at the end of this file also record the actual memory use of the
process (including threads), as reported by the operating system at the end of
the step, and the total memory still in use per label. Note this includes
memory still active from previous steps and the total memory is also continued
from the previous dump.

MPI task communication reports
------------------------------

When SWIFT is configured using the ``--enable-mpiuse-reports`` flag it will
log any all asynchronous MPI communications made to send particle updates
between nodes to support the tasks.

The output files are called ``mpiuse_report-rank<m>-step<n>.dat``, i.e. one
per rank per step. These have a line for each request for communication, either
an MPI_Irecv or MPI_Isend and a line for the subsequent completion (successful
MPI_Test).

Each line of the logs contains the following information:

.. code-block:: none

   stic:             ticks since the start of this step
   etic:             ticks since the start of the simulation
   dtic:             ticks that the request was active
   step:             current step
   rank:             current rank
   otherrank:        rank that the request was sent to or expected from
   type itype:       task type as string and enum
   subtype isubtype: task subtype as string and enum
   activation:       1 if record for the start of a request, 0 if request completion
   tag:              MPI tag of the request
   size:             size, in bytes, of the request
   sum:              sum, in bytes, of all requests that are currently not logged as complete 

The stic values should be synchronized between ranks as all ranks have a
barrier in place to make sure they start the step together, so should be
suitable for matching between ranks. The unique keys to associate records
between ranks (so that the MPI_Isend and MPI_Irecv pairs can be identified)
are "otherrank/rank/subtype/tag/size" and "rank/otherrank/subtype/tag/size"
for send and recv respectively. When matching ignore step0.




Task and Threadpool Plots and Analysis Tools
--------------------------------------------

A variety of plotting tools for tasks and threadpools is available in ``tools/task_plots/``.
To be able to use the task analysis tools, you need to compile swift with ``--enable-task-debugging``
and then run swift with ``-y <interval>``, where ``<interval>`` is the interval between time steps
on which the additional task data will be dumped. Swift will then create ``thread_stats-step<nr>.dat``
and ``thread_info-step<nr>.dat`` files. Similarly, for threadpool debugging, you need to compile
swift with ``--enable-threadpool-debugging`` and then run it with ``-Y <interval>``.

For the analysis and plotting scripts listed below, you need to provide the **\*info-step<nr>.dat** 
files as a cmdline argument, not the ``*stats-step<nr>.dat`` files.

A short summary of the scripts in ``tools/task_plots/``:

- ``analyse_tasks.py``: 
    The output is an analysis of the task timings, including deadtime per thread
    and step, total amount of time spent for each task type, for the whole step
    and per thread and the minimum and maximum times spent per task type.
- ``analyse_threadpool_tasks.py``: 
    The output is an analysis of the threadpool task timings, including 
    deadtime per thread and step, total amount of time spent for each task type, for the
    whole step and per thread and the minimum and maximum times spent per task type.
- ``iplot_tasks.py``: 
    An interactive task plot, showing what thread was doing what task and for 
    how long for a step.  **Needs python2 and the tkinter module**.
- ``plot_tasks.py``: 
    Creates a task plot image, showing what thread was doing what task and for how long. 
- ``plot_threadpool.py``: 
    Creates a threadpool plot image, showing what thread was doing what threadpool call and for
    how long. 


For more details on the scripts as well as further options, look at the documentation at the top
of the individual scripts and call them with the ``-h`` flag.
