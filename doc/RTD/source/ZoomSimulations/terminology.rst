.. Zoom Terminology
   Will Roper, 20th November 2024

Terminology
===========

Zoom Region
~~~~~~~~~~~

The zoom region is a cubic region encompassing the high-resolution particles of a zoom-in resimulation.

Zoom cell
~~~~~~~~~

Zoom cells are high-resolution cells tesselating the zoom region. Zoom cells can get any task type but can only have ``regular`` cell subtypes. 

Background cell
~~~~~~~~~~~~~~~

Background cells are low-resolution cells which tesselate the entire volume (including the zoom region). Background cells only ever get gravity (and related common) tasks. They can have ``void`` or ``neighbour`` cell subtypes.

Buffer cell
~~~~~~~~~~~

Buffer cells are an (optional) intermediate resolution cell type which pad the volume between the background and zoom cells. Buffer cells are only used when a zoom region is sufficiently small relative to the parent box to need a buffer around the zoom region. Including buffer avoids pathological numbers of background cells or poorly resolving the zoom region. As with background cells, buffer cells only get gravity (and related common) tasks and can have 


Void cell
~~~~~~~~~

Background (and buffer cells if used) must align with the bounds of the zoom region. These background/buffer cells above the zoom region are denoted “void” cells, they will never be popualted by particles explicitly. They do however have multipoles for long-range and multipole-multipole gravity interactions.

Neighbour cell 
~~~~~~~~~~~~~~ 

A Neighbour cell can be either a background or a buffer cell (if they are included). These are cells within the maximum distance at which we will have to do a non-mesh gravity interaction (e.g. a long-range gravity task, ``grav_mm`` or ``pair/grav`` task).
