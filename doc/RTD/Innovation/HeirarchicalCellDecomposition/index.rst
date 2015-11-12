.. _GettingStarted:

Heirarchical Cell Decomposition
=================================

Most SPH codes rely on spatial trees to decompose the simulation space. This decomposition makes neighbour-finding simple, at the cost of computational efficiency. Neighbour-finding using the tree-based approach has an average computational cost of ~O(logN) and has a worst case behaviour of ~O(N^2/3), both cases grow with the total number of particles N. SWIFT's neighbour-finding algorithm however, has a constant scaling of ~O(1) per particle. This results from the way SWIFT decomposes its domain.

The space is divided up into a grid of rectangular cells with an edge length that is greater than or equal to the maximum smoothing of any particle in the simulation, h_max (See Fig. 1). 

.. figure:: InitialDecomp.pdf
   :height: 100px
   :width: 200 px
   :scale: 50 %
   :alt: 2D Cell Decomposition
   :align: center

In this initial decomposition if a particle p_j is within range of particle p_i, both will either be in the same cell (self-interaction) or in neighbouring cells (pair-interaction). Each cell then only has to compute its self-interactions and pair-interactions for each of its particles.   

The best case scenario is when each cell only contains particles that have a smoothing length equal to the cell edge length and even then, for any given particle p_i it will only interact with 16% of the total number of particles in the same cell and surrounding neighbours. This percentage decreases if the cell contains particles whose smoothing length is less than the cell edge length. Therefore the cell decomposition needs to be refined recursively by bisecting a cell along each dimension if the following conditions are met:

1) The cell contains more than a minimum number of particles

2) The smoothing length of a reasonable number of particles within a cell is less than half the cell's edge length

Once the cells have been split their self-interactions can be decomposed into self-interactions of its sub-cells and corresponding pair interactions (See Fig. 2).

.. figure:: SplitCell.pdf
   :height: 100px
   :width: 200 px
   :scale: 50 %
   :alt: Refined Cell Decomposition
   :align: center


.. toctree::
   :maxdepth: 1
