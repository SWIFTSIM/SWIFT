.. External potentials in SWIFT
   Folkert Nobels, 25th October 2018
   
External Potentials 
===================

SWIFT can be run with an external potential on this page we will summarize the
current potentials which can be run with SWIFT and how to implement your own 
potential in SWIFT.

Implemented External Potentials
-------------------------------

Currently there are several potentials implemented in SWIFT. On this page we 
give a short overview of the potentials that are implemented in the code:

1. No potential (none)
2. Point mass potential (point-mass): classical point mass, can be placed at
   a position with a mass.
3. Plummer potential (point-mass-softened): in the code a softened point mass 
   corresponds to a Plummer potential, can be placed at a position with a mass.
4. Isothermal potential (isothermal): An isothermal potential which corresponds 
   to a density profile which is :math:`\propto r^{-2}` and a potential which is 
   logarithmic. This potential has as free parameters the rotation velocity 
   and the position.
5. Hernquist potential (hernquist): A potential that is given by the Hernquist 
   potential: 
   
   :math:`\Phi(r) = - \frac{GM}{r+a}.`

   The free parameters of Hernquist potential are mass, scale length,
   and softening. The potential can be set at any position in the box.
6. NFW potential (nfw): The most used potential to describe dark matter halos, the  
   potential is given by:

   :math:`\Phi(r) = - \frac{4\pi G \rho_0 R_s^3}{r} \ln \left( 1+ 
   \frac{r}{R_s} \right).`

   This potential has as free parameters the concentration of the DM halo, the
   virial mass (:math:`M_{200}`) and the critical density.
7. Sine wave (sine-wave)
8. Point mass ring (point-mass-ring)
9. Disc Patch (disc-patch)


How to implement your own potential
-----------------------------------

The first step in implementing your own potential is making a directory of your
potential in the ``src/potential`` folder and creating a file in the folder 
called ``potential.h``.

Configuring the potential 
^^^^^^^^^^^^^^^^^^^^^^^^^

To get started you can copy a ``potential.h`` file from an already implemented 
potential. In this potential the header guards (e.g. ``#IFDEF <>``) need to be 
changed to the specific potential and the ``struct`` and 
``potential_init_backend`` need to be  changed such that it uses your potential 
and reads the correct potential from the parameter file during running the 
program.

Add the potential to the ``potential.h`` file in the ``src`` directory such that
the program knows that it is possible to run with this potential.

Furthermore during the configuration of the code it also needs to be clear for 
the program that the code can be configured to run with the different 
potentials. This means that the ``configure.ac`` file needs to be changed.
This can be done to add an other case in the potential::

  case "$with_potential" in
     none)
        AC_DEFINE([EXTERNAL_POTENTIAL_NONE], [1], [No external potential])
     ;;
     newpotential)
        AC_DEFINE([EXTERNAL_POTENTIAL_NEWPOTENTIAL], [1], [New external potential])
     ;;

After this change it is possible to configure the code to use your new potential.

