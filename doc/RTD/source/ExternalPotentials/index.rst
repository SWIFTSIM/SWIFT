.. External potentials in SWIFT
   Folkert Nobels, 25th October 2018
   Alejandro Benitez-Llambay, October 2019
   Matthieu Schaller, December 2021
   
External Potentials 
===================

SWIFT can be run with an external potential on this page we will summarize the
current potentials which can be run with SWIFT and how to implement your own 
potential in SWIFT.

Implemented External Potentials
-------------------------------

Currently there are several potentials implemented in SWIFT. On this page we
give a short overview of the potentials that are implemented in the code. They
are all switched on at configuration time via the argument
``--with-ext-potential=XYZ`` and get their own independant section in the
parameter file. The name of the potential to pass to the configuration script is
given in the parenthesis of the titles below; for instance
``--with-ext-potential=herquist``.

1. No potential (``none``)
^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the default setup. No external potential is used.
There are no parameters associated with this model.

2. Constant acceleration (``constant``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This "potential" just applies a constant acceleration vector at every position
:math:`\vec{x}` in the simulation volume. This is very common in idealised test
cases like the Rayleigh-Taylor instability or engineering applications relying
on Earth's constant gravity field.

 * :math:`\phi(\vec{x}) = -\vec{g} \cdot \vec{x}`
 * :math:`\vec{a}(\vec{x}) = \vec{g}`

The only parameter of this model is the vector :math:`\vec{g}` given in `cgs`
units:

.. code:: YAML
	  
   ConstantPotential:
      g_cgs:           [0., 0., -9,81]  # Earth acceleration along z-axis (cgs units)

      
3. Point mass potential (``point-mass``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simple single point mass that can be placed at any position :math:`\vec{p}`.
The position can be specified either absolutely or with respect to the centre of
the simulation volume. 

 * :math:`\phi(\vec{x}) = -\displaystyle\frac{G_{\rm N}m}{|r|}`
 * :math:`\vec{a}(\vec{x}) = -\displaystyle\frac{G_{\rm N}m}{|r|^3}\vec{r}`,

where :math:`\vec{r} = \vec{x} - \vec{p}`.

The code also imposes an extra limit of the size of the particles'
time-step. The time-step has to be shorter than :math:`\Delta t_{pot} = C
|\vec{a}(\vec{x})| / |\dot{\vec{a}}(\vec{x})|`. This criterion is designed to
make sure the changes in accelerations remain small. The jerk
(:math:`\dot{\vec{a}}\equiv\frac{d}{dt}\vec{a}`) is calculated from the
positions of the particles and only includes the contribution of the external
potential itself. The other criteria (CFL, self-gravity, ...) are applied on top
of this criterion. The dimensionless constant :math:`C` defaults to `FLT_MAX` if
left unspecified, effectively making SWIFT run without this time-step criterion.

The parameters of the model are:

.. code:: YAML
	  
   PointMassPotential:
      position:         [3., 4., 5.]  # Location of the point mass (internal units)
      useabspos:                   1  # Use absolute positions (0 for relative to centre)
      mass:                 5.972e24  # Mass of the point (internal units)
      timestep_mult:             0.1  # (Optional) The dimensionless constant C in the time-step condition
	  
4. Plummer potential (``point-mass-softened``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A single point mass with a fixed softening length using a Plummer shape that can
be placed at any position :math:`\vec{p}`.  The position can be specified either
absolutely or with respect to the centre of the simulation volume.

 * :math:`\phi(\vec{x}) = -\displaystyle\frac{G_{\rm N}m}{\sqrt{|r|^2 + \epsilon^2}}`
 * :math:`\vec{a}(\vec{x}) = -\displaystyle\frac{G_{\rm N}m}{(|r|^2 + \epsilon^2)^{3/2}}\vec{r}`,

where :math:`\vec{r} = \vec{x} - \vec{p}`.

The code also imposes an extra limit of the size of the particles'
time-step. The time-step has to be shorter than :math:`\Delta t_{pot} = C
|\vec{a}(\vec{x})| / |\dot{\vec{a}}(\vec{x})|`. This criterion is designed to
make sure the changes in accelerations remain small. The jerk
(:math:`\dot{\vec{a}}\equiv\frac{d}{dt}\vec{a}`) is calculated from the
positions of the particles and only includes the contribution of the external
potential itself. The other criteria (CFL, self-gravity, ...) are applied on top
of this criterion. The dimensionless constant :math:`C` defaults to `FLT_MAX` if
left unspecified, effectively making SWIFT run without this time-step criterion.

The parameters of the model are:

.. code:: YAML
	  
   PointMassPotential:
      position:         [3., 4., 5.]  # Location of the point mass (internal units)
      useabspos:                   1  # Use absolute positions (0 for relative to centre)
      mass:                 5.972e24  # Mass of the point (internal units)
      softening:                0.01  # Softening length (internal units)
      timestep_mult:             0.1  # (Optional) The dimensionless constant C in the time-step condition
      
5. Isothermal potential (``isothermal``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An isothermal potential which corresponds to a density profile which is
:math:`\propto r^{-2}` and a potential which is logarithmic. This potential is
entirely set by its centre :math:`\vec{p}` and the (radially-constant) rotation
velocity of the potential :math:`v_{\rm rot}`. The centre of the potential is softened.
	 
 * :math:`\phi(\vec{x}) = -\displaystyle\frac{1}{4\pi G_{\rm N}}\log(\sqrt{|r|^2 + \epsilon^2})`
 * :math:`\vec{a}(\vec{x}) = -\displaystyle\frac{v_{\rm rot}^2} {|r|^2 + \epsilon^2}`,

where :math:`\vec{r} = \vec{x} - \vec{p}`.

The code also imposes an extra limit of the size of the particles'
time-step. The time-step has to be shorter than :math:`\Delta t_{pot} = C
|\vec{a}(\vec{x})| / |\dot{\vec{a}}(\vec{x})|`. This criterion is designed to
make sure the changes in accelerations remain small. The jerk
(:math:`\dot{\vec{a}}\equiv\frac{d}{dt}\vec{a}`) is calculated from the
positions of the particles and only includes the contribution of the external
potential itself. The other criteria (CFL, self-gravity, ...) are applied on top
of this criterion. The dimensionless constant :math:`C` defaults to `FLT_MAX` if
left unspecified, effectively making SWIFT run without this time-step criterion.

The parameters of the model are:

.. code:: YAML
	  
   IsothermalPotential:
      position:         [3., 4., 5.]  # Location of the centre of the profile (internal units)
      useabspos:                   1  # Use absolute positions (0 for relative to centre)
      vrot:                      200  # Rotation velocity of the profile (internal units)
      softening:                0.01  # Softening length (internal units)
      timestep_mult:             0.1  # (Optional) The dimensionless constant C in the time-step condition

   
6. Hernquist potential (``hernquist``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can set up a potential as given by Hernquist (1990): 

:math:`\Phi(r) = - \frac{G_{\rm N}M}{r+a},`

where :math:`M` is the total Hernquist mass and :math: `a` is the Hernquist-
equivalent scale radius of the potential. The potential can be set at any 
position in the box. It adds an additional time step constraint that limits
the time step to a maximum of a specified fraction of the circular orbital 
time at the current radius of the particle. The other criteria 
(CFL, self-gravity, ...) are applied on top of this criterion. For example, a
fraction of 0.01 means that an orbital period would be resolved by 100 time steps.

In the most basic version, the Hernquist potential can be run by providing 
only the Hernquist mass, scale length, softening length and fraction of the 
orbital time for the time stepping. In this case the model parameters may 
look something like:

.. code:: YAML

    HernquistPotential:
        useabspos:       0              # 0 -> positions based on centre, 1 -> absolute positions 
        position:        [0.,0.,0.]     # Location of centre of isothermal potential with respect to centre of the box (if 0) otherwise absolute (if 1) (internal units)
        mass:            1e12           # Mass of the Hernquist potential
        scalelength:     2.0            # scale length a (internal units)
        timestep_mult:   0.01           # Dimensionless pre-factor for the time-step condition, determines the fraction of the orbital time we use to do the time integration; fraction of 0.01 means we resolve an orbit with 100 timesteps
        epsilon:         0.2            # Softening size (internal units)

Besides the basic version, it is also possible to run the Hernquist 
potential for idealised disk galaxies that follow the approach of 
`Springel, Di Matteo & Hernquist (2005)
<https://ui.adsabs.harvard.edu/abs/2005MNRAS.361..776S/abstract>`_. This 
potential, however, uses a corrected value of the formulation that improves 
the match with the NFW profile (below) with the same M200 (Nobels+ in prep). 
Contrary to the above (idealizeddisk: 0 setup), the idealised disk setup runs 
by specifying one out of :math:`M_{200}`, :math:`V_{200}`, or :math:`R_{200}`, 
plus concentration and reduced Hubble constant.

In this case, we don't provide the 'mass' and 'scalelength' parameters, but
'M200' (or 'V200', or 'R200') and 'concentration' :math:`c`, as well as reduced Hubble
constant :math:`h` to define the potential. The parameters of the model may look something like:

.. code:: YAML

    HernquistPotential:
        useabspos:       0              # 0 -> positions based on centre, 1 -> absolute positions 
        position:        [0.,0.,0.]     # Location of centre of isothermal potential with respect to centre of the box (if 0) otherwise absolute (if 1) (internal units)
        idealizeddisk:   1              # Run with an idealized galaxy disk
        M200:            137.0          # M200 of the galaxy disk
        h:               0.704          # reduced Hubble constant (value does not specify the used units!)
        concentration:   9.0            # concentration of the Halo
        diskfraction:    0.040          # Disk mass fraction
        bulgefraction:   0.0            # Bulge mass fraction
        timestep_mult:   0.01           # Dimensionless pre-factor for the time-step condition, determines the fraction of the orbital time we use to do the time integration; fraction of 0.01 means we resolve an orbit with 100 timesteps
        epsilon:         0.2            # Softening size (internal units)

The user should specify one out of 'M200', 'V200', or 'R200' to define
the potential. The reduced Hubble constant is then used to determine the
other two. Then, the scale radius is calculated as :math:`R_s = R_{200}/c`,
where :math:`c` is the concentration, and the Hernquist-equivalent scale-length
is calculated as:

:math:`a = \frac{b+\sqrt{b}}{1-b} R_{200}`,

where :math:`b = \frac{2}{c^2}(\ln(1+c) - \frac{c}{1+c})`.

Two examples using the Hernquist potential can be found in ``swiftsim/examples/GravityTests/``. 
The ``Hernquist_radialinfall`` example puts 5 particles with zero velocity in a Hernquist 
potential, resulting in radial orbits. The ``Hernquist_circularorbit``example puts three
particles on a circular orbit in a Hernquist potential, one at the inner region, one at the
maximal circular velocity, and one in the outer region. To run these examples, SWIFT must
be configured with the flag ``--with-ext-potential=hernquist``, or ``hernquist-sdmh05``
(see below).

7. Hernquist SDMH05 potential (``hernquist-sdmh05``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the same potential as Hernquist with the difference being the
way that the idealised disk part is calculated. This potential uses
exactly the same approach as `Springel, Di Matteo & Hernquist (2005)
<https://ui.adsabs.harvard.edu/abs/2005MNRAS.361..776S/abstract>`_,
which means that ICs generated with the original `MakeNewDisk` code can
be used with this potential. Contrary to the updated Hernquist
potential (above), it is not possible to have an identically matched
NFW potential in this case.

The parameters needed for this potential are the same set of variables as 
above, i.e. 'mass' and 'scalelength' when we don't use the idealised
disk, and 'concentration' plus one out of 'M200', 'V200', or 'R200' if 
we do. As one of the three is provided, the reduced Hubble constant is
used to calculate the other two. Then, the scale radius is calculated
using the NFW definition, :math:`R_s = R_{200}/c`, and the Hernquist-
equivalent scale length is given by

:math:`a = R_s \sqrt{2(\ln(1+c) - \frac{c}{1+c})}.`

Runs that provide e.g. M200 and c (using idealised disk) are thus equivalent
to providing mass and scale length if calculated by the above equation (without
idealized disk). 

   
8. Navarro-Frenk-White potential (``nfw``):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
The most used potential to describe dark matter halos, the  
potential is given by:

:math:`\Phi(r) = - \frac{4\pi G_{\rm N} \rho_0 R_s^3}{r} \ln \left( 1+ 
\frac{r}{R_s} \right).`

This potential has as free parameters the concentration of the DM halo, the
virial mass (:math:`M_{200}`) and the critical density. The potential add 
an additional time step constrain that limits the time step to a maximum of 
a specified fraction of the orbital time. 

This potential in the centre and the enclosed mass at R200 are identical to 
the Hernquist potential discussed above. 

The parameters of the model are:

.. code:: YAML

    NFWPotential:
        useabspos:       0              # 0 -> positions based on centre, 1 -> absolute positions 
        position:        [0.,0.,0.]     # Location of centre of isothermal potential with respect to centre of the box (if 0) otherwise absolute (if 1) (internal units)
        M200:            137.0          # M200 of the galaxy disk
        h:               0.704          # reduced Hubble constant (value does not specify the used units!)
        concentration:   9.0            # concentration of the Halo
        diskfraction:    0.040          # Disk mass fraction
        bulgefraction:   0.0            # Bulge mass fraction
        timestep_mult:   0.01           # Dimensionless pre-factor for the time-step condition, basically determines the fraction of the orbital time we use to do the time integration
        epsilon:         0.2            # Softening size (internal units)

   
9. NFW poential + Miyamoto-Nagai potential (``nfw-mn``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   This includes and NFW potential (identical to nfw)
   plus an axisymmetric Miyamoto-Nagai potential. The Miyamoto-Nagai potential is given by:

   :math:`\Phi(R,z) = - \frac{G_{\rm N} M_{d}}{\sqrt{R^2 + \left ( R_d + \sqrt{z^2 + Z_d^2} \right )^2}}`,

   where :math:`R^2 = x^2 + y^2` is the projected radius and :math:`M_d`, :math:`R_d`, :math:`Z_d` are the 
   mass, scalelength and scaleheight of the disk (in internal units), respectively. 

The parameters of the model are:

.. code:: YAML

    NFW_MNPotential:
        useabspos:        0              # 0 -> positions based on centre, 1 -> absolute positions 
        position:         [0.,0.,0.]     # Location of centre of isothermal potential with respect to centre of the box (if 0) otherwise absolute (if 1) (internal units)
        M_200:            137.0          # M200 of the halo in internal units
	critical_density: 123.4          # Critical density of the universe in internal units
        concentration:    9.0            # concentration of the NFW halo
	Mdisk:            3.0            # Mass of the disk in internal units
	Rdisk:            3.0            # Disk scale-length in internal units
	Zdisk:            3.0            # Disk scale-height in internal units
        timestep_mult:    0.01           # Dimensionless pre-factor for the time-step condition, basically determines the fraction of the orbital time we use to do the time integration

   
10. Sine wave (``sine-wave``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This "potential" is designed for specific idealised tests. It is a basic (not
really physical!) sinusoid wave along the x-axis with a unit wavelength and a
tunable amplitude and growth time.

 * :math:`\phi(\vec{x}) = A \cos\left(2 \pi x_{0}\right) / 2\pi`
 * :math:`a_0(\vec{x}) = A \sin\left(2 \pi x_{0}\right)`,

where the 0 subscript indicates the x-component. The y- and z-components are zero.

The amplitude :math:`A` can be growing linearly to its maximum over a fixed
time-scale.

Optionally, a constant maximail time-step size can be used with this potential.

The parameters of the model are:

.. code:: YAML

   SineWavePotential:
      amplitude:                 2.5  # The amplitude of the wave (internal units) 
      growth_time:               1.2  # Time for the amplitude to grow to its final value (internal units)
      timestep_limit:            5.0  # (Optional) The max time-step of the particles (internal units)

11. Disc Patch (``disc-patch``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A potential corresponding to a vertical slice through a galactic disk. This
follows the definitions of `Creasey, Theuns & Bower (2013)
<https://adsabs.harvard.edu/abs/2013MNRAS.429.1922C>`_ equations (16) and (17).
The potential is implemented along the x-axis.

12. MWPotential2014 (``MWPotential2014``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This potential is based on ``galpy``'s ``MWPotential2014`` from `Jo Bovy (2015) <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B>`_ and consists in a NFW potential for the halo, an axisymmetric Miyamoto-Nagai potential for the disk and a bulge modeled by a power spherical law with exponential cut-off. The bulge is given by the density:

:math:`\rho(r) = A \left( \frac{r_1}{r} \right)^\alpha \exp \left( - \frac{r^2}{r_c^2} \right)`,

where :math:`A` is an amplitude, :math:`r_1` is a reference radius for amplitude, :math:`\alpha` is the inner power and :math:`r_c` is the cut-off radius.

The resulting potential is:

:math:`\Phi_{\mathrm{MW}}(R, z) = f_1 \Phi_{\mathrm{NFW}} + f_2 \Phi_{\mathrm{MN}} + f_3 \Phi_{\text{bulge}}`,

where :math:`R^2 = x^2 + y^2` is the projected radius and :math:`f_1`, :math:`f_2` and :math:`f_3` are three coefficients that adjust the strength of each individual component.

The parameters of the model are:

.. code:: YAML

    MWPotential2014Potential:
      useabspos:        0          # 0 -> positions based on centre, 1 -> absolute positions
      position:         [0.,0.,0.] # Location of centre of potential with respect to centre of the box (if 0) otherwise absolute (if 1) (internal units)
      timestep_mult:    0.005      # Dimensionless pre-factor for the time-step condition, basically determines the fraction of the orbital time we use to do the time integration
      epsilon:          0.001      # Softening size (internal units)
      concentration:    9.823403437774843      # concentration of the Halo
      M_200_Msun:       147.41031542774076e10  # M200 of the galaxy disk (in M_sun)
      H:                1.2778254614201471     # Hubble constant in units of km/s/Mpc
      Mdisk_Msun:       6.8e10                 # Mass of the disk (in M_sun)
      Rdisk_kpc:        3.0                    # Effective radius of the disk (in kpc)
      Zdisk_kpc:        0.280                  # Scale-height of the disk (in kpc)
      amplitude:        1.0                    # Amplitude of the bulge
      r_1_kpc:          1.0                    # Reference radius for amplitude of the bulge (in kpc)
      alpha:            1.8                    # Exponent of the power law of the bulge
      r_c_kpc:          1.9                    # Cut-off radius of the bulge (in kpc)
      potential_factors: [0.4367419745056084, 1.002641971008805, 0.022264787598364262] #Coefficients that adjust the strength of the halo (1st component), the disk (2nd component) and the bulge (3rd component)

Note that the default value of the "Hubble constant" here seems odd. As it
enters multiplicatively with the :math:`f_1` term, the absolute normalisation is
actually not important.
      
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

