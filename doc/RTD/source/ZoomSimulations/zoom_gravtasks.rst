.. Zoom Gravity Task Creation
   Will Roper, 25th September 2024

Gravity Task Creation
=====================

To reduce the number of pair tasks between cell grids and in general reduce the number of zoom->zoom pair tasks, we need to offload as many interactions as possible into multipole interactions. To achieve this we first create tasks at the highest possible level and then split these tasks checking if we can do multipole interactions as we split.

Top Level Task Creation
-----------------------

In `zoom_maketasks` we create the initial gravity tasks for each background cell just as would be done in a uniform box. Background cells get self tasks and pair tasks. This will include void cells (the cells that contain the zoom region) which get self and pair tasks despite containing no particles. Tasks involving void cells will be split as described below.

If buffer cells are being used then the void cells are instead buffer cells and these also get self and pair tasks. In addition to the buffer->buffer pair tasks (which will include the void cell tasks) we also explictly create buffer->background pair tasks to handle the interactions between these cell grids since they are not nested in the same way as the zoom cells.

Void Task Splitting
-------------------

Once all tasks are created at the top level (and between background and buffer cells if the latter are being used) we begin the process of recursively splitting the tasks down the cell trees. This is true for all self or pair tasks to exceed the number of particles designated for splitting in the parameter file. For void cells this logic is different though. We have created top level tasks in void cells that don't actually contain the particles, we need to split there tasks until we reach at least the top level of the zoom cells (which contain the particles).

We always split any tasks involving void cells until we reach a level without void cells. As we split, we can check to see if the interaction can instead be done via multipoles, in which case a ``grav_mm`` task is made. If a ``grav_mm`` task cannot be used (because the cells are too close together) then we recurse down to the next level in the void cell tree/s. Once we hit the zoom cell leaves the task splitting calls the normal logic based on task size is used.

For a void->void interaction we can freely split both sides of the interaction because we know there will be progeny (void cells are by construction going to always have 8 progeny). However, if there is only 1 void cell in the interaction we first split the void side of the interaction until we hit a non-void cell. Once we have two non-void cells we can split as normal.

This process guarantees we have the most chances to avoid a direct interaction and interact at the highest possible level across the zoom region boundary.

Void Tasks
----------

Void cells can only have a limited subset of tasks. Since they can't contain particles they can't have self and pair tasks (hence the explicit splitting). Void cells at the super level can have:

- ``grav_init``: The multipoles need to be initialised for each step they are active.
- ``grav_long_range``: Long range interactions need to be done at the same level the gravity tasks are done. If there is no void super level (i.e. no ``grav_mm`` tasks) then the zoom cell progeny will instead get long range interactions.
- ``grav_mm``: Cheap multipole interactions can be done in the void tree instead of numerous direct pair interactions at the zoom level.
- ``grav_down``: Since the multipoles can be interacted the potentials need to be propagated down to the zoom cells. The void cell ``grav_down`` will handle the propagation down to the zoom super level. The zoom ``grav_down`` will then handle the propagation from there.

The drifting of void cells is handled outside of the tasking prior to the call to ``engine_launch``.
