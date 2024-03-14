.. Zoom Hydrodynamics
   Will Roper, 14th March 2024

Hydrodynamics with a Zoom Region
================================

Baryons, and therefore the hydrodynamics tasks, are isolated to the high resolution zoom region. This means baryons that leave the zoom region will no longer have their properties updated.

This is another of the many reasons why padding the zoom region beyond the extent of the high resolution particles is so important. Choose too small a value for `ZoomRegion:region_pad_factor` and baryons will leave the zoom region regularly. Choose too large a value and the zoom region will be poorly resolved losing the efficiency gains from introducing the higher resolution zoom cells.

Pressureless Boundary
---------------------

Only having high resolution baryons introduces a pressureless boundary at the edge of the high resolution particle distribution...


Maximum Smoothing Length
~~~~~~~~~~~~~~~~~~~~~~~~
