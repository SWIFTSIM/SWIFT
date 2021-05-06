.. Neutrinos
   Willem Elbers, 7 April 2021

.. _neutrinos:

Neutrino implementation
=======================

SWIFT can also accurately model the effects of massive neutrinos in
cosmological simulations. At the background level, massive neutrinos
and other relativistic species can be included by specifying their
number and masses in the cosmology section of the parameter file
(see :ref:`Parameters_cosmology`).

At the perturbation level, neutrinos can be included as a separate particle
species (``PartType6``). To facilitate this, SWIFT implements the
:math:`\delta f` method for shot noise suppression (`Elbers et al. 2020
<https://ui.adsabs.harvard.edu/abs/2020arXiv201007321E/>`_). The method
works by statistically weighting the particles during the simulation,
with weights computed from the Liouville equation using current and
initial momenta. The method can be activated by specifying
``Neutrino:use_delta_f`` in the parameter file.

The implementation of the :math:`\delta f` method in SWIFT assumes a
specific method for generating the initial neutrino momenta (see below).
This makes it possible to reproduce the initial momentum when it is
needed without increasing the memory footprint of the neutrino particles.
If perturbed initial conditions are not needed, the initial momenta can
be generated internally by specifying ``Neutrino:generate_ics`` in the
parameter file. This will assign ``PartType6`` particles to each
neutrino mass specified in the cosmology and generate new velocities
based on the homogeneous (unperturbed) Fermi-Dirac distribution.

Relativistic Drift
------------------

At high redshift, neutrino particles move faster than the speed of light
if the usual Newtonian expressions are used. To rectify this, SWIFT
implements a relativistic drift correction. In this convention, the
internal velocity variable (see theory/Cosmology) is
:math:`v^i=a^2u^i=a^2\dot{x}^i\gamma^{-1}`, where :math:`u^i` is the
spatial part of the 4-velocity, :math:`a` the scale factor, and
:math:`x^i` a comoving position vector. The conversion factor to the
coordinate 3-velocity is :math:`\gamma=ac/\sqrt{a^2c^2+v^2}`. This
factor is applied to the neutrino particles throughout the simulation.

Generating Fermi-Dirac momenta
------------------------------

The implementation of the :math:`\delta f` method in SWIFT assumes that
neutrinos were initially assigned a Fermi-Dirac momentum using the following
method. Each particle has a fixed 64-bit unsigned integer :math:`\ell` given
by the particle ID [#f1]_ (plus an optional seed: ``Neutrino:neutrino_seed``).
This number is transformed into a floating point number :math:`u\in(0,1)`,
using the following pseudo-code based on splitmix64:

.. code-block:: none

    m = l + 0x9E3779B97f4A7C15
    m = (m ^ (m >> 30)) * 0xBF58476D1CE4E5B9;
    m = (m ^ (m >> 27)) * 0x94D049BB133111EB;
    m = m ^ (m >> 31);
    u = (m + 0.5) / (UINT64_MAX + 1)

This is subsequently transformed into a Fermi-Dirac momentum
:math:`q = F^{-1}(u)` by evaluating the quantile function. To generate
neutrino particle initial conditions with perturbations, one first generates
momenta from the unperturbed Fermi-Dirac distribution using the above method
and then applies perturbations in any suitable manner.

When using the :math:`\delta f` method, SWIFT also assumes that ``PartType6``
particles are assigned to all :math:`N_\nu` massive species present in the
cosmology, such that the particle with fixed integer :math:`\ell` corresponds
to species :math:`i = \ell\; \% \;N_\nu\in[0,N_\nu-1]`.

The sampled Fermi-Dirac speeds and neutrino masses are written into the
snapshot files as ``SampledSpeeds`` and ``MicroscopicMasses``.

.. [#f1] Currently, it is not guaranteed that a particle ID is unique.
