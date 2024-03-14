.. Zoom Terminology
   Will Roper, 14th March 2024

Terminology
===========

Zoom Region
-----------

The zoom region is a cubic region encompassing the high-resolution particles of a zoom-in resimulation.

Zoom cell
---------

Zoom cells are high-resolution cells tesselating the zoom region. Zoom cells can get any task type.

Background cell
---------------

Background cells are low-resolution cells which tesselate the entire volume (including the zoom region). Background cells only ever get gravity (and related common) tasks.

Buffer cell
-----------

Buffer cells are an intermediate resolution cell type which pad the volume between the background and zoom cells. Buffer cells are only used when a zoom region is sufficiently small relative to the parent box to need a buffer around the zoom region so as to avoid pathological numbers of background cells. As with background cells, buffer cells only get gravity (and related common) tasks.


Void cell
---------
Background/buffer cells must align with the bounds of the zoom region with 1 or more background/buffer cells in the zoom region. These background/buffer cells in the zoom region are “void” cells, in the sense that they will never get tasks or explicit particles. They do however have multipoles for long-range gravity calculations.

Empty cell
----------

An empty cell is a background cell above the buffer region (around the zoom region). Empty cells are only used when buffer cells are also used. An empty cell can never contain anything or get any kind of task. They exist purely to maintain the cell array layout.
