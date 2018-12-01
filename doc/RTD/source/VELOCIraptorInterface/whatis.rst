.. What is VELOCIraptor
   Folkert Nobels 12th October 2018


What is VELOCIraptor?
=====================

.. toctree::    
   :maxdepth: 2    
   :hidden:    
   :caption: Contents: 

In SWIFT it is possible to run a cosmological simulation and at the same time
do on the fly halo finding at specific predefined intervals. For finding the 
halos SWIFT uses VELOCIraptor (Elahi, Thacker and Widrow; 2011) [#velo]_, this 
is a C++ halo finder that can use MPI. It differs from other halo finder 
algorithms in the sense that it uses the velocity distributions of the 
particles in the simulations and the the positions of the particles to get
a better estimate of which particles are part of a specific halo and 
whether there are substructures in halos. 

The Algorithm
-------------

The VELOCIraptor algorithm consists basically of the following steps [#ref]_:

1. A kd-tree is constructed based on the maximization of the Shannon-entropy,
   this means that every level in the kd-tree an equal number of particles 
   are distributed between the 8 lower nodes. This is based on their position
   and their corresponding density, this results in more equal density 
   distributed nodes. This is also the implicit step in the algorithm that 
   takes into account the absolute positions of the particles.
2. The next part is calculating the the centre of mass velocity and the 
   velocity distribution for every individual node in the kd-tree. 
3. Then the algorithm estimates the background velocity density function for
   every particle based on the cell of the particle and the six nearest
   neighbour cells. This prevents the background velocity density function 
   to be over sensitive for variations between different cells due to dominant
   halo features in the velocity density function. 
4. After this the algorithm searches for the nearest velocity neighbours 
   (:math:`N_v`) from a set of nearest position neighbours (:math:`N_x>N_v`).
   The neighbours' positions do not need to be in the cell of the particles, in
   general the set of nearest position neighbours is substantially larger than
   the nearest velocity neighbours, the default is set as :math:`N_x=32 N_v`.
5. The individual local velocity density function is calculated for every 
   particle.
6. The fractional difference is calculated between the local velocity density 
   function and the background velocity density function.
7. Based on the calculated ratio, outliers are picked and the outliers are  
   grouped together in halos and subhalos.
  


.. Every halo finder has limitations, the limitations of VELOCIraptor are:

.. 1. The algorithm is mostly sensitive to substructures that are on the tail
   of the Gaussian velocity density function, this means that VELOCIraptor
   is most sensitive for subhalos which are cold (slow rotating) but have 
   a large bulk velocity


.. _Velociraptor: http://adsabs.harvard.edu/abs/2011MNRAS.418..320E
.. [#velo] For technical information regarding VELOCIraptor see: Velociraptor_
.. [#ref] This part is based on the explanation given in the Elahi, Thacker and
   Widrow (2011) paper (Velociraptor_)
