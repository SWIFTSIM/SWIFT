.. AnalysisTools
   Loic Hausammann 20th March 2019
   Peter W. Draper 28th March 2019
   Mladen Ivkovic 18th March 2021
   Bert Vandenbroucke 31st February 2022

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


Task dependencies for a single cell
-----------------------------------

There is an option to additionally write the dependency graphs of the task dependencies for a single cell.
You can select which cell to write using the ``Scheduler:dependency_graph_cell: cellID`` parameter, where ``cellID`` is the cell ID of type long long.
This feature will create an individual file for each step specified by the ``Scheduler:dependency_graph_frequency`` and, differently from the full task graph, will create an individual file for each MPI rank that has this cell.

Using this feature has several requirements:

- You need to compile SWIFT including either ``--enable-debugging-checks`` or ``--enable-cell-graph``. Otherwise, cells won't have IDs.
- There is a limit on how many cell IDs SWIFT can handle while enforcing them to be reproducibly unique. That limit is up to 32 top level cells in any dimension, and up to 16 levels of depth. If any of these thresholds are exceeded, the cells will still have unique cell IDs, but the actual IDs will most likely vary between any two runs.

To plot the task dependencies, you can use the same script as before: ``tools/plot_task_dependencies.py``. The dependency graph now may have some tasks with a pink-ish background colour: These tasks represent dependencies that are unlocked by some other task which is executed for the requested cell, but the cell itself doesn't have an (active) task of that type itself in that given step.


Task levels
-----------------

At the beginning of each simulation the file ``task_level_0.txt`` is generated.
It contains the counts of all tasks at all levels (depths) in the tree.
The depths and counts of the tasks can be plotted with the script ``tools/plot_task_levels.py``.
It will display the individual tasks at the x-axis, the number of each task at a given level on the y-axis, and the level is shown as the colour of the plotted point.
Additionally, the script can write out in brackets next to each task's name on the x-axis on how many different levels the task exists using the ``--count`` flag.
Finally, in some cases the counts for different levels of a task may be very close to each other and overlap on the plot, making them barely visible.
This can be alleviated by using the ``--displace`` flag:
It will displace the plot points w.r.t. the y-axis in an attempt to make them better visible, however the counts won't be exact in that case.
If you wish to have more task level plots, you can use the parameter ``Scheduler:task_level_output_frequency``.
It defines how many steps are done in between two task level output dumps.




Cell graph
----------

An interactive graph of the cells is available with the configuration option ``--enable-cell-graph``. During a
run, SWIFT will generate a ``cell_hierarchy_*.csv`` file per MPI rank at the frequency given by the parameter
``--cell-dumps=n``. The script ``tools/make_cell_hierarchy.py`` can be used to collate the files produced by
different MPI ranks and convert them into a web page that shows an interactive cell hierarchy. The script
takes the names of all the files you want to include as input, and requires an output prefix that will be used
to name the output files ``prefix.csv`` and ``prefix.html``. If the prefix path contains directories that do
not exist, the script will create those.

The output files cannot be directly viewed from a browser, because they require a server connection to
interactively load the data. You can either copy them over to a server, or set up a local server yourself. The
latter can also be done directly by the script by using the optional parameter ``--serve``.

When running a large simulation, the data loading may take a while (a few seconds for EAGLE_6). Your browser
should not be hanging, but will appear to be idle. For really large simulations, the browser will give up and
will probably display an error message.

If you wish to add some information to the graph, you can do it by modifying the files ``src/space.c`` and
``tools/data/cell_hierarchy.html``. In the first one, you will need to modify the calls to ``fprintf`` in the
functions ``space_write_cell_hierarchy`` and ``space_write_cell``. Here the code is simply writing CSV files
containing all the required information about the cells. In the second file, you will need to find the
function ``mouseover`` and add the field that you have created. You can also increase the size of the bubble
through the style parameter ``height``.

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

The stic values should be synchronised between ranks as all ranks have a
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
and ``thread_info-step<nr>.dat`` files. Similarly, for threadpool related tools, you need to compile
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

Task data is also dumped when using MPI and the tasks above can be used on
that as well, some offer the ability to process all ranks, and others to
select individual ranks.

It is also possible to process a complete run of task data from all the
available steps using the ``process_plot_tasks.py`` and
``process_plot_tasks_MPI.py`` scripts, as appropriate.
These scripts have one required argument: a time limit to use on the horizontal
time axis. When set to 0, this limit is determined by the data for each step,
making it very hard to compare relative sizes of different steps.
The optional ``--files`` arguments allows more control over which steps are
included in the analysis. Large numbers of tasks can be analysed more
efficiently by using multiple processes (the optional ``--nproc`` argument),
and if sufficient memory is available, the parallel analysis can be optimised
by using the size of the task data files to schedule parallel processes more
effectively (the ``--weights`` argument).


.. _dumperThread:

Live internal inspection using the dumper thread
------------------------------------------------

If the configuration option ``--enable-dumper`` is used then an extra thread
is created that polls for the existence of local files called
``.dump<.rank>``. When found this will trigger dump logs of the current state
of various internal queues and loggers, depending on what is enabled.

Without any other options this will dump logs of the current tasks in the
queues (these are those ready to run when time and all conflicts allow) and
all the tasks that are expected to run this step (those which are active in
the current time step). If ``memuse-reports`` is enabled the currently logged
memory use is also dumped and if ``mpiuse-reports`` is enabled the MPI
communications performed this step are dumped. As part of this dump a report
about MPI messages which have been logged but not completed is also made to
the terminal. These are useful when diagnosing MPI deadlocks.

The active tasks are dumped to files ``task_dump-step<n>.dat`` or
``task_dump_MPI-step<n>.dat_<rank>`` when using MPI.

Similarly the currently queued tasks are dumped to files
``queue_dump-step<n>.dat`` or ``queue_dump_MPI-step<n>.dat_<rank>``.

Memory use logs are written to files ``memuse-error-report-rank<n>.txt``.
The MPI logs follow the pattern using ``mpiuse-error-report-rank<n>.txt``.

The ``.dump<.rank>`` files once seen are deleted, so dumping can be done more
than once. For a non-MPI run the file is simply called ``.dump``, note for MPI
you need to create one file per rank, so ``.dump.0``, ``.dump.1`` and so on.


Deadlock Detector
---------------------------

When configured with ``--enable-debugging-checks``, the parameter

.. code-block:: yaml

    Scheduler:
        deadlock_waiting_time_s:   300.

can be specified. It specifies the time (in seconds) the scheduler should wait
for a new task to be executed during a simulation step (specifically: during a
call to ``engine_launch()``). After this time passes without any new tasks being
run, the scheduler assumes that the code has deadlocked. It then dumps the same
diagnostic data as :ref:`the dumper thread <dumperThread>` (active tasks, queued
tasks, and memuse/MPIuse reports, if swift was configured with the corresponding
flags) and aborts.

A value of zero or a negative value for ``deadlock_waiting_time_s`` disable the
deadlock detector.

You are likely well advised to try and err on the upper side for the time to
choose for the ``deadlock_waiting_time_s`` parameter. A value in the order of
several (tens of) minutes is recommended. A too small value might cause your run to
erroneously crash and burn despite not really being deadlocked, just slow or
badly balanced.






Neighbour search statistics
---------------------------

One of the core algorithms in SWIFT is an iterative neighbour search
whereby we try to find an appropriate radius around a particle's
position so that the weighted sum over neighbouring particles within
that radius is equal to some target value. The most obvious example of
this iterative neighbour search is the SPH density loop, but various
sub-grid models employ a very similar iterative neighbour search. The
computational cost of this iterative search is significantly affected by
the number of iterations that is required, and it can therefore be
useful to analyse the progression of the iterative scheme in detail.

When configured with ``--enable-ghost-statistics=X``, SWIFT will be
compiled with additional diagnostics that statistically track the number
of iterations required to find a converged answer. Here, ``X`` is a
fixed number of bins to use to collect the required statistics
(``ghost`` refers to the fact that the iterations take place inside the
ghost tasks). In practice, this means that every cell in the SWIFT tree
will be equipped with an additional ``struct`` containing three sets of
``X`` bins (one set for each iterative neighbour loop: hydro, stellar
feedback, AGN feedback). For each bin ``i``, we store the number of
particles that required updating during iteration ``i``, the number of
particles that could not find a single neighbouring particle, the
minimum and maximum smoothing length of all particles that required
updating, and the sum of all their search radii and all their search
radii squared. This allows us to calculate the upper and lower limits,
as well as the mean and standard deviation on the search radius for each
iteration and for each cell. Note that there could be more iterations
required than the number of bins ``X``; in this case the additional
iterations will be accumulated in the final bin. At the end of each time
step, a text file is produced (one per MPI rank) that contains the
information for all cells that had any relevant activity. This text file
is named ``ghost_stats_ssss_rrrr.txt``, where ``ssss`` is the step
counter for that time step and ``rrrr`` is the MPI rank.

The script ``tools/plot_ghost_stats.py`` takes one or multiple
``ghost_stats.txt`` files and computes global statistics for all the
cells in those files. The script also takes the name of an output file
where it will save those statistics as a set of plots, and an optional
label that will be displayed as the title of the plots. Note that there
are no restrictions on the number of input files or how they relate;
different files could represent different MPI ranks, but also different
time steps or even different simulations (which would make little
sense). It is up to the user to make sure that the input is actually
relevant.
