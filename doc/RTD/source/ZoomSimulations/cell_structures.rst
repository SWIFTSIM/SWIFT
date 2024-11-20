.. Zoom cell structures
   Will Roper, 14th March 2024

Cell Construction and Their Hierarchy
=====================================

At it's root SWIFT is a cartesian grid of cells holding particles, on which computations can be performed. When defining a single cartesian grid for a zoom simulation the calculation becomes inbalanced with most of the work confined to a single (or handful of) cell(s). To combat this a hierachy of cell grids is introduced to better resolve the work. 

There are two possible heirarchies that can be used depending on the geometry of the zoom simulation: (i) a two part hierachy with background cells and zoom cells nested at a certain depth inside the background cells, (iii) a three part hierachy with background cells, buffer cells, and zoom cells. Buffer cells are nested at a certain depth inside the background cells, with zoom cells then nested at a further depth inside the buffer cells. The former is used when the zoom region is large relative to the parent box, the latter when the zoom region is small relative to the parent box. 

Ultimately which method is used is automatically chosen by SWIFT based on the parameters defined the zoom region geometry and the high resolution particle distribution. If the zoom region can be aligned with the edges of background cells without adding too much padding around the high-resolution particles then only background and zoom cells are needed. Buffers cells will automatically be added when the padding to align background and zoom cells increases the zoom region size by a factor of 2. 

Large Zoom Regions
------------------

.. figure:: figures/zoom_geometry_nobuffer.png
            :width: 400px
            :align: center
            :alt: Large zoom region cells


Small Zoom Regions
------------------

.. figure:: figures/small_cells.png
            :width: 400px
            :align: center
            :alt: Small zoom region cells


Cell Construction
-----------------

In any of the methods detailed above the process of cell construction is:

1. Define zoom region geometry and background cell properties.
2. Shift particles to centre the zoom region.
3. Construct zoom, background, and (when applicable) buffer top level cell grids.
4. Recursively construct cell trees and multipoles in all top level cells.
5. (In MPI land) Communicate multipoles.
6. Construct the void cell trees out of the void cells and zoom multipoles.


Below is a general description of the important stages in the construction of the cell hierarchy.

Zoom Region Dimensions
~~~~~~~~~~~~~~~~~~~~~~

The dimensions of the zoom region are defined in ``zoom_init.zoom_region_init`` by first finding the extent of all non-background particles and then padding this extent by ``ZoomRegion:region_pad_factor`` (defined in the parameter file, by default 1.5). After this initial definition, the zoom region dimensions can be increased to ensure the background cell grid/s align (based on each method detailed above).

The number of zoom cells along an axis of the zoom region is defined by the ``ZoomRegion:zoom_top_level_cells`` parameter, while their size is calculated in ``zoom_init.zoom_region_init`` after the dimensions of the zoom region have been found. If running with hydro, only zoom cells are given hydro-related tasks.

Particle Shifting
~~~~~~~~~~~~~~~~~

Before constructing cells we shift the particle distribution to place the geometric midpoint of the high-resolution particle distribution at the centre of the box. This is done to ensure boundary effects can be ignored while constructing the cell grids and tasks.

This shift is independant of the user specified shift defined in the parameter file (``InitialConditions:shift``). The shift applied to centre the zoom region will be undone prior to writing out any positions to a snapshot. This is not true of ``InitialConditions:shift`` which will be respected and not undone.

Void Cell Tree
~~~~~~~~~~~~~~

Once the 2 (or 3 when buffer cells are included) cell grids in the hierachy have been constructed, a cell tree is constructed in all void cells. This cell tree is constructed in the same way as a “normal” cell tree is for any other cell. However, a void cell tree’s leaves are not leaves of the void cell but are instead zoom cells (which themselves contain a normal cell tree).

To avoid recursing from a void cell right through to the leaves of the zoom cell tree the parent void cells of the zoom cells are given ``c->split = 0``, while zoom cell's parents are ``NULL`` but have a new member (``void_parent``) which points to their parent void cell in the void cell tree.

This linking of zoom cells as leaves is one of the reasons the cell grids must align perfectly (although task definitions and proxies also rely on this being the case). Since the leaves must be linked in like this ``ZoomRegion:zoom_top_level_cells`` must be ``ZoomRegion:region_buffer_cell_ratio`` times a power of two.

The void cell trees allow for long-range gravity tasks involving the zoom region to be done at levels above individual zoom cells and thus limits the number of long-range gravity calculations around the zoom region. They also provide a method for limiting the number of MPI communications around the zoom region.
