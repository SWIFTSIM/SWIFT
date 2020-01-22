.. GEAR sub-grid model
   Matthieu Schaller, 20th December 2018


GEAR model
===========

Pressure Floor
~~~~~~~~~~~~~~

In order to avoid the artificial collapse of unresolved clumps, a minimum in pressure is applied to the particles.
This additional pressure can be seen as the pressure due to unresolved hydrodynamics turbulence and is given by:

.. math::
    P_\textrm{Jeans} = \frac{\rho}{\gamma} \left( \frac{4}{\pi} G h^2 \rho N_\textrm{Jeans}^{2/3} - \sigma^2 \right)

where :math:`\rho` is the density, :math:`\gamma` the adiabatic index, :math:`G` is the gravitational constant,
:math:`h` the smoothing length, :math:`N_\textrm{Jeans}` is the number of particle required in order to resolve a clump and
:math:`\sigma` the velocity dispersion.


This must be directly implemented into the hydro schemes, therefore only a subset of schemes (Gadget-2, SPHENIX and Pressure-Energy) have the floor available.
In order to implement it, you need equation 12 in `Hopkins 2013 <https://arxiv.org/abs/1206.5006>`_:

.. math::
   m_i \frac{\mathrm{d}v_i}{\mathrm{d}t} = - \sum_j x_i x_j \left[ \frac{P_i}{y_i^2} f_{ij} \nabla_i W_{ij}(h_i) + \frac{P_j}{y_j^2} f_{ji} \nabla_j W_{ji}(h_j) \right]

and simply replace the :math:`P_i, P_j` by the pressure with the floor.
Here the :math:`x, y` are simple weights that should never have the pressure floor included even if they are related to the pressure (e.g. pressure-entropy).


Cooling: Grackle
~~~~~~~~~~~~~~~~
   
Grackle is a chemistry and cooling library presented in `B. Smith et al. 2016 <https://arxiv.org/abs/1610.09591>`_ 
(do not forget to cite if used).  Four different modes are available:
equilibrium, 6 species network (H, H\\( ^+ \\), e\\( ^- \\), He, He\\( ^+ \\)
and He\\( ^{++} \\)), 9 species network (adds H\\(^-\\), H\\(_2\\) and
H\\(_2^+\\)) and 12 species (adds D, D\\(^+\\) and HD).  Following the same
order, the swift cooling options are ``grackle_0``, ``grackle_1``, ``grackle_2``
and ``grackle_3`` (the numbers correspond to the value of
``primordial_chemistry`` in Grackle).  It also includes some self-shielding
methods and UV background.  In order to use the Grackle cooling, you will need
to provide an HDF5 table computed by Cloudy.

When starting a simulation without providing the different fractions, the code
supposes an equilibrium and computes the fractions automatically.

In order to compile SWIFT with Grackle, you need to provide the options ``with-chemistry=GEAR`` and ``with-grackle=$GRACKLE_ROOT``
where ``$GRACKLE_ROOT`` is the root of the install directory (not the ``lib``).

You will need a Grackle version later than 3.1. To compile it, run
the following commands from the root directory of Grackle:
``./configure; cd src/clib``.
Update the variables ``LOCAL_HDF5_INSTALL`` and ``MACH_INSTALL_PREFIX`` in
the file ``src/clib/Make.mach.linux-gnu``.
Finish with ``make machine-linux-gnu; make && make install``.
If you encounter any problem, you can look at the `Grackle documentation <https://grackle.readthedocs.io/en/latest/>`_

You can now provide the path given for ``MACH_INSTALL_PREFIX`` to ``with-grackle``.
