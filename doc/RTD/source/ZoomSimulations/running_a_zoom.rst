.. Running a zoom
   Will Roper, 14th March 2024

Running a Zoom Simulation
=========================

Very little extra is needed to run a zoom simulation. All you need is the command line option, a set of zoom initial conditions, and the ``ZoomRegion`` parameter file block. Each of these is described below.

Command Line Option
-------------------

To tell SWIFT to run with a zoom region you simply include ``--zoom`` (or ``-z``) at runtime. This option can be combined with almost all other options, i.e. to run a zoom with the eagle model and 128 threads you can invoke ``swift --zoom --eagle -t 128 <path/to/params.yml>``.


Initial Conditions
------------------

The initial conditions should have all hydrodynamic particles stored under ``PartType0`` and all high resolution dark matter particles should be stored under ``PartType1``. All background dark matter particles should be stored under ``PartType2``. For further details on the initial conditions see their :doc:`docs <../InitialConditions/index>`. No extra data is needed in the initial conditions to run a zoom beyond the background particles.


Parameter File
--------------

The ``ZoomRegion`` parameter block is used to define the properties of the zoom region and background. An example is shown below.

.. code-block:: yaml

  ZoomRegion:
    region_pad_factor:          1.1
    bkg_top_level_cells:        16
    zoom_top_level_depth:       2
    buffer_top_level_depth:     0
    neighbour_max_tree_depth:   2
    bkg_subdepth_diff_grav:     3

The top 4 parameters here control the construction of the cell hierarchy enabling the efficient division of a zoom calculation (see the :doc:`cell docs <cell_structures>` for more information), while the latter 2 control the construction of tasks within the cell tree.

Below is a description of each of these parameters.

``region_pad_factor``
~~~~~~~~~~~~~~~~~~~~~

 The extent of the zoom region is defined as the extent of the high resolution particle distribution mutliplied by ``region_pad_factor``. This is not to be confused with any padding applied during initial condition generation, that is naturally included in the high-resolution particle distribution. The ``region_pad_factor`` should be selected to be as small as possible while still providing a sufficient padding when running with hydrodynamics (see :doc:`Zoom Hydrodynamics <zoom_hydrodynamics`). **If the padding is set too large then the zoom region will be poorly resolved and performance will be sub-optimal.**

 It's important to note that the value of ``region_pad_factor`` eventually used during the simulation will be modified slightly during cell construction to ensure the zoom region tesselates the parent volume correctly. The padding factor after tesselation is reported by the code if run with verbosity enabled (``-v 1`` at runtime). Should the requested cell geometry lead to a substantial increase in the padding factor then an error will be thrown, at which point you should review the parameters you have set. You probably need to increase the background resolution and/or modify the cell depths (more below).

``bkg_top_level_cells``
~~~~~~~~~~~~~~~~~~~~~~~~

This replaces ``Scheduler:max_top_level_cells`` (which is ignored when running a zoom simulation). It defines the number of background cells along each axis of the whole box.

This value is always respected during cell construction (although it can be later modified in a regrid if the maximum smoothing length in the zoom region demands it).


``zoom_top_level_depth``
~~~~~~~~~~~~~~~~~~~~~~~~ 

This defines the depth of the zoom cells in a single background cell (inside the cell tree). For instance, a zoom cell depth of 4 would mean a single background cell will contain ``2**4 = 16`` zoom cells at the leaves of its cell tree.

This value is always respected during cell construction (although it can be later modified in a regrid if the maximum smoothing length in the zoom region demands it and increasing the number of background cells isn't sufficient).

``buffer_top_level_depth`` 
~~~~~~~~~~~~~~~~~~~~~~~~~~

This parameter is only required when the zoom region is small relative background cells and buffer cells are needed to pad the space inbetween the background cell boundaries and the zoom region. This defines the depth of the buffer cells in a single background cell (inside the cell tree). For instance, a buffer cell depth of 2 would mean a single background cell will contain ``2**2 = 4`` buffer cells at the leaves of its cell tree. By default it is 0 and the code will error and tell you if you need to define a non-zero value and haven't. ``buffer_top_level_depth`` must always be less than ``zoom_top_level_depth``.

This value is always respected during cell construction and is not modified during a regrid.

``neighbour_max_tree_depth``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

This parameter controls the maximum depth neighbour cells will be **forced** to split to. Splitting neighbour cells deeper can enable better division of interactions across the zoom region boundary. Note that this only controls how far neighbour cells will be forced to be split to if they weren't going to be split anyway, if a neighbour cell naturally needs splitting beyond this depth because of the number of particles it contains that will still be split.  

``bkg_subdepth_diff_grav`` 
~~~~~~~~~~~~~~~~~~~~~~~~~~ 

This defines the difference between the depth of the cell leaves and the level at which we will define gravity tasks for background and buffer cells. In a normal simulation this is controlled by ``cell_subdepth_diff_grav``, but in a zoom simulation its useful to differentiate this parameter for background and zoom cells. The normal ``cell_subdepth_diff_grav`` parameter is used for zoom cells while this one covers all other cell types.
