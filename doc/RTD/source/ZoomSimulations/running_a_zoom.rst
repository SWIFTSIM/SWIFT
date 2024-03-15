.. Running a zoom
   Will Roper, 14th March 2024

Running a Zoom Simulation
=========================

To run a zoom simulation all you need is the command line option, a set of zoom initial conditions, and the ``ZoomRegion`` parameter file block. Each of these is described below.

Command Line Option
-------------------

To tell SWIFT to run with a zoom region you simply include ``--zoom`` (or ``-z``) at runtime. This option can be combined with almost all other options, i.e. to run a zoom with the eagle model and 128 threads you can invoke ``swift --zoom --eagle -t 128 <path/to/params.yml>``.


Initial Conditions
------------------

The initial conditions should have all hydrodynamic particles stored under ``PartType0`` and all high resolution dark matter particles should be stored under ``PartType1``. All background dark matter particles should be stored under ``PartType2``. For further details on the initial conditions see their :doc:`docs <../InitialConditions/index>`. No other extra information is needed in the initial conditions to run a zoom.


Parameter File
--------------

The ``ZoomRegion`` parameter block is used to define the properties of the zoom region and background. An example is shown below.

.. code-block:: none

  ZoomRegion:
    region_pad_factor: 1.1
    zoom_top_level_cells: 8
    bkg_top_level_cells: 8
    region_buffer_cell_ratio: 2

These 4 parameters alone control the construction of the cell hierarchy which enables efficient zoom simulations in SWIFT (see the :doc:`cell docs <cell_structures>` for more information).

``region_pad_factor``
~~~~~~~~~~~~~~~~~~~~~

 The extent of the zoom region is defined as the extent of the high resolution particle distribution mutliplied by ``region_pad_factor``. This is independant of any padding applied during initial condition generation and should be selected to be as small as possible while still providing a sufficient padding when running with hydrodynamics (see :doc:`Zoom Hydrodynamics <zoom_hydrodynamics`). If the padding is set too large then the zoom region will be poorly resolved and performance will be sub-optimal.

 It's important to note that the value of ``region_pad_factor`` will be modified slightly during cell construction to ensure the zoom region tesselates the parent volume correctly in each cell scenario (:doc:`discussed here <cell_structures>`). The padding factor after tesselation is reported by the code. Should the requested cell geometry lead to a substantial increase in the padding factor then an error will be thrown, at which point you should review the parameters you have set. You probably need to increase the background resolution and/or include buffer cells (more below).

``zoom_top_level_cells``
~~~~~~~~~~~~~~~~~~~~~~~~

This replaces ``Scheduler:max_top_level_cells`` (which is ignored when running a zoom) and defines the number of zoom cells along each axis of the zoom region.

This value is always respected during cell construction (although it can be later modified in a regrid if the hydrodynamics demand it).

``bkg_top_level_cells``
~~~~~~~~~~~~~~~~~~~~~~~

This is the background cell analogue of ``Scheduler:max_top_level_cells`` and defines the number of background cells along an axis of the full parent volume.

This value behaves differently in different scenarios:

- For large zoom regions: This value is treated as a target to aim for.
- For intermediate zoom regions: This value is ignored entirely and the number of background cells is the number of zoom regions which tesselate the parent volume along each axis.
- For small zoom regions: This value is respected and never modified.


``region_buffer_cell_ratio``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This defines the number of buffer cells that tesselate the zoom region and should only be used for small zoom regions. When set to 0 buffer cells won't be made (i.e. there will only be background and zoom cells).

When ``region_buffer_cell_ratio > 0`` ``zoom_top_level_cells`` should be a power of 2 multiplied by ``region_buffer_cell_ratio`` to ensure the cell edges align. If this is not the case an error will be thrown. If ``region_buffer_cell_ratio`` is non-zero in a situation it shouldn't be an error will also be thrown.
