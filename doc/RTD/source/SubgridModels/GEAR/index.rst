.. GEAR sub-grid model
   Matthieu Schaller, 20th December 2018


GEAR model
===========


Cooling: Grackle
~~~~~~~~~~~~~~~~
   
Grackle is a chemistry and cooling library presented in `B. Smith et al. 2016 <https://arxiv.org/abs/1610.09591>`_ 
(do not forget to cite if used).  Four different modes are available:
equilibrium, 6 species network (H, H\\( ^+ \\), e\\( ^- \\), He, He\\( ^+ \\)
and He\\( ^{++} \\)), 9 species network (adds H\\(^-\\), H\\(_2\\) and
H\\(_2^+\\)) and 12 species (adds D, D\\(^+\\) and HD).  Following the same
order, the swift cooling options are ``grackle``, ``grackle1``, ``grackle2``
and ``grackle3`` (the numbers correspond to the value of
``primordial_chemistry`` in Grackle).  It also includes some self-shielding
methods and UV background.  In order to use the Grackle cooling, you will need
to provide an HDF5 table computed by Cloudy.

When starting a simulation without providing the different fractions, the code
supposes an equilibrium and computes the fractions automatically.

In order to compile SWIFT with Grackle, you need to provide the options ``with-grackle``
and ``with-chemistry``.

You will need a Grackle version later than 3.1. To compile it, run
the following commands from the root directory of Grackle:
``./configure; cd src/clib``.
Update the variables ``LOCAL_HDF5_INSTALL`` and ``MACH_INSTALL_PREFIX`` in
the file ``src/clib/Make.mach.linux-gnu``.
Finish with ``make machine-linux-gnu; make && make install``.
If you encounter any problem, you can look at the `Grackle documentation <https://grackle.readthedocs.io/en/latest/>`_

You can now provide the path given for ``MACH_INSTALL_PREFIX`` to ``with-grackle``.
