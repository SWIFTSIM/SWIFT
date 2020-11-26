.. Planetary Initial Conditions
   Jacob Kegerreis, 13th March 2020

.. _planetary_initial_conditions:
   
Initial Conditions
==================

Tools for the creation of initial conditions are available via 
the 
`WoMa <https://github.com/srbonilla/WoMa>`_ and 
`SEAGen <https://github.com/jkeger/seagen>`_ open-source python packages, 
including: creating spherical or spinning planetary (or similar) profiles;
placing particles to match arbitrary profiles with precise SPH densities;
and setting the initial target and impactor positions and velocities,
as presented in 
`Kegerreis et al. (2019)  <https://doi.org/10.1093/mnras/stz1606>`_ and
Ruiz-Bonilla et al. (2020).

They are available with documentation and examples at 
https://github.com/srbonilla/WoMa and https://github.com/jkeger/seagen,
or can be installed directly with ``pip``
(https://pypi.org/project/woma/, https://pypi.org/project/seagen/).


Settling initial conditions with fixed entropies
------------------------------------------------

If the particles' equations of state include specific entropies, 
and the initial conditions file includes specific entropies for each particle
(in ``PartType0/Entropies``), 
then configuring SWIFT with ``--enable-planetary-fixed-entropy``
will override the internal energy of each particle each step such that its 
specific entropy remains constant. 

This should be used with caution, but may be a convenient way to maintain an 
entropy profile while initial conditions settle to equilibrium with their 
slightly different SPH densities.
