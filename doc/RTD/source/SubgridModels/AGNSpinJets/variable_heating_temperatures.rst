.. AGN spin and jet model
   Filip Husko, 26 September 2023

.. AGN_spin_jet:

Variable AGN heating temperature model
--------------------------------------

As part of this AGN model we also introduce an option of using variable AGN heating temperatures, as opposed to a constant value. This can be turned on by setting ``AGN_heating_temperature_model`` to ``AGN_heating_temperature_local``. In this scheme, we use three different criteria, all of which are applied simultaneously, and all of which are numerically motivated. These are: 1) the `Dalla Vecchia & Schaye (2012) <https://ui.adsabs.harvard.edu/abs/2012MNRAS.426..140D/abstract>`_ to prevent numerical overcooling, 2) a replenishment condition to make sure that the BH smoothing kernel can be replenished at a rate similar to the rate of evacuation due to heating and 3) a crossing condition that ensures a sufficient separation of heating events in time such that the kernel is never expected to contain more than a single heated particle.

We implement the `Dalla Vecchia & Schaye (2012) <https://ui.adsabs.harvard.edu/abs/2012MNRAS.426..140D/abstract>`_ condition using Eqn. (18) from the paper. We invert it to obtain the minimum heating temperature:

.. math::
    \Delta T_\mathrm{DV}= 1.8\times10^6\hspace{0.5mm}\mathrm{K}\hspace{0.5mm}\times\bigg(\frac{N_\mathrm{ngb}m_\mathrm{g}}{60\times10^6\hspace{0.5mm}\mathrm{M}_\odot} \bigg)^{1/3}\bigg(\frac{n_\mathrm{H}}{0.1\hspace{0.5mm}\mathrm{cm}^{-3}} \bigg)^{2/3},

where :math:`N_\mathrm{ngb}` is the number of neighbours in the BH kernel, :math:`m_\mathrm{g}` the average neighbour mass and :math:`n_\mathrm{H}` the number density of gas in the BH kernel.

We derive the second condition by comparing the time-scales of replenishment onto the kernel and evacuation out of it on account of gas heating. The evacuation time-scale is the time required to evacuate the entire kernel of gas, and it is given by :math:`t_\mathrm{evac}=m_\mathrm{g}N_\mathrm{ngb}/\dot{M}`, where :math:`\dot{M}` is the mass flux associated with feedback. It can be expressed using the feedback power :math:`P` and internal energy per unit mass :math:`e` as :math:`\dot{M}=P/e`, where :math:`e=3k_\mathrm{B}\Delta T/2\mu m_\mathrm{p}`, :math:`\Delta T` is the heating temperature, :math:`m_\mathrm{p}` the proton mass and :math:`\mu\approx0.6` the mean molecular weight for ionized gas.

The replenishment time-scale is given by :math:`t_\mathrm{repl}=h/v_\mathrm{repl}`, where :math:`v_\mathrm{repl}` is an effective velocity with which gas can flow inwards to replenish the kernel. Given this formulation, we can set the two time-scales (of evacuation and replenishment) equal to each other and solve for a heating temperature that ensures timely replenishment:

.. math::
    \Delta T_\mathrm{repl} = \frac{2\mu m_\mathrm{p}}{3 k_\mathrm{B}}\frac{hP}{N_\mathrm{ngb}m_\mathrm{g}v_\mathrm{repl}}.

Finally, we calculate the replenishment velocity :math:`v_\mathrm{repl}` as follows. We assume that gas can replenish the kernel under the effects of either gas turbulence or pressure. We thus express the replenishment velocity as :math:`v_\mathrm{repl} = \max(\sigma,\Tilde{c}_\mathrm{s})`, where :math:`\sigma` is the velocity dispersion of gas and :math:`\Tilde{c}_\mathrm{s}` an effective sound speed, which we calculate by assuming that there is always some ISM exherting its pressure on gas near the BH, even if all of the gas in the kernel is cold. We choose the form :math:`\Tilde{c}_\mathrm{s}=\max(c_\mathrm{s,hot},\hspace{0.5mm}10\hspace{0.5mm}\mathrm{km}\mathrm{s}^{-1})`, where :math:`c_\mathrm{s,hot}` is the kernel-weighted average sound speed of all particles that have a temperature :math:`T>10^4` K and :math:`10` km :math:`\mathrm{s}^{-1}` is the sound speed of the ISM assuming a typical temperature of :math:`T=10^4` K.

The third and final condition we use is based on the time it takes a single heated particle to cross and exit the kernel before the next one is heated. The time-interval between two heating events is :math:`\Delta = m_\mathrm{g}/\dot{M}`, while the time required for a heated particle to cross the kernel is :math:`\Delta t_\mathrm{cross}= h/c_\mathrm{s,\Delta T}`, where :math:`c_\mathrm{s,\Delta T} = \sqrt{\gamma(\gamma-1)e} = \sqrt{5k_\mathrm{B}\Delta T/3\mu m_\mathrm{p}}` is the sound speed of the heated gas. Equating these two time-scales, we obtain the final heating temperature:

.. math::
    \Delta T_\mathrm{cross} = \frac{\mu m_\mathrm{p}}{k_\mathrm{B}}\bigg(\frac{2hP}{\sqrt{15}m_\mathrm{g}}\bigg)^{2/3}.

The final heating temperature scheme we use can be written as:

.. math::
    \Delta T = \max[\xi \max(\Delta T_\mathrm{DV}, \Delta T_\mathrm{repl}, \Delta T_\mathrm{cross}), \Delta T_\mathrm{min}],

where :math:`\Delta T_\mathrm{min}` is an additional temperature floor meant to prevent very-low temperature heating events, and :math:`\xi` is an additional free parameter that can be used to rescale heating temperatures from all three conditions to higher values, which may be necessary to correctly calibrate the simulations. Specifically, one may use the hot gas fractions to choose a value of :math:`\xi` if :math:`\Delta T_\mathrm{min}` is set to a low value of :math:`\approx10^{7.5}` K, or one may set :math:`\xi=1`, thus using the numerical conditions exactly as derived, and instead calibrate the simulations by varying :math:`\Delta T_\mathrm{min}`.

