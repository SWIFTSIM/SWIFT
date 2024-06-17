.. AGN spin and jet model
   Filip Husko, 1 April 2022

.. AGN_spin_jet:

Jet launching algorithm
-----------------------

In order to launch jets, we introduce a jet reservoir that functions identically to the thermal reservoir used in EAGLE and COLIBRE. When the jet reservoir exceeds the value :math:`N_\mathrm{j}\overline{m}_\mathrm{ngb}v_\mathrm{j}^2/2`, kicks are given out to :math:`N_\mathrm{j}` particles. Here, :math:`N_\mathrm{j}` is the target number of particles to kick per each kicking event (typically equal to 2, but should always be a multiple of 2 to ensure approximately symmetrical jets). :math:`\overline{m}_\mathrm{ngb}v_\mathrm{j}^2/2` is the average kinetic energy of one kicking event, with :math:`\overline{m}_\mathrm{ngb}` the mean neighbour mass of gas particles in the BH smoothing kernel and :math:`v_\mathrm{j}` the target kicking velocity.

These kicks are handed out to particles in a symmetric way with respect to the spin vector of the BH. :math:`N_\mathrm{j}/2` particles are kicked from the 'upper' hemisphere relative to the spin vector, and the other half from the lower hemisphere. The particles to be kicked can be any in the smoothing kernel. We include four different choices: the particles kicked are: 1) the closest to the BH, 2) the farthest from the BH, 3) the ones of minimal density and 4) the ones closest to the spin axis, in terms of angular distance. Note that these sortings are done for each hemisphere seperately. 

The particles chosen are always given velocities based on the same algorithm, regardless of their positions in the kernel. We perform the actual kicks in the following way. Velocity kicks are chosen at random from a cone around the current spin vector with a (half-)opening angle of :math:`\theta_\mathrm{j}`. In particular, we first choose the kick vector around the z-axis as :math:`\hat{v}_\mathrm{kick}=(\sin\theta\cos\phi,\hspace{0.3mm}\sin\theta\sin\phi,\hspace{0.3mm}\cos \theta)`. Here, :math:`\cos\theta` is chosen uniformly from the interval :math:`[\cos\theta_\mathrm{j},1]`, and :math:`\sin\theta=\sqrt{1-\cos\theta^2}`. :math:`\phi` is chosen uniformly from :math:`[0,2\pi]`. This random vector, now representing a random kick within a cone around the z-axis, is rotated into the frame of the spin vector so that the cone is pointing in the right direction. For particles being kicked from the 'negative' side of the BH hemisphere, the final kick vector is simply multiplied by :math:`-1`.

We increase the particle's velocity in the chosen kick direction by an amount :math:`\vert \Delta \vec{v} \vert` such that its energy increases by :math:`(1/2)m_\mathrm{gas}v_\mathrm{jet}^2`. For this reason, :math:`\vert \Delta \vec{v} \vert< v_\mathrm{jet}` generally holds. We calculate :math:`\vert \Delta \vec{v} \vert` from the equation

.. math::
    (\vec{v}_0+\Delta\vec{v})^2 = \vert \vec{v}_0\vert ^2 + v_\mathrm{j}^2,
    
which follows from conservation of energy. This vector equation can be solved to yield the necessary magnitude of the velocity increase that should be applied (in the chosen kick direction)

.. math::
    \vert \Delta\vec{v}\vert = \sqrt{v_\mathrm{j}^2 + v_\mathrm{0,j}^2} - v_\mathrm{0,j},
    
where :math:`v_\mathrm{0,j} = \sin \theta\vert \vec{v}_0\vert` is the magnitude of the initial velocity projected onto the chosen kick direction, with :math:`\sin \theta` the angle between the direction of the initial velocity and the chosen kick direction.

Black hole time steps
---------------------

Given the changes to the BH model in the form of AGN jets and BH spin evolution, a few additional time-step criteria need to be implemented. The minimum of these time-steps is taken to actually evolve the BH, alongside the other time-steps already used for the BH in the code. We introduce a jet-related time-step that is given by:

.. math::
    \Delta t_\mathrm{jet}=\frac{\Delta E_\mathrm{jet}}{P_\mathrm{jet}}.

This time-step ensures that the BH is woken up by the time it needs to 'hand out' a pair of kicks. In the above equation, :math:`P_\mathrm{jet}` is the current, instantenous jet power, while :math:`\Delta E_\mathrm{jet}=2\times m_\mathrm{ngb}v_\mathrm{jet}^2` is the energy to be handed out to a pair of particles, with :math:`m_\mathrm{ngb}` the average gas particle mass in the BH's kernel, and :math:`v_\mathrm{jet}` the target jet velocity.

We also introduce two time-steps related to the angular momentum of the BH. The first of these ensures that the magnitude of spin does not change too much over a single time-step, and it is given by

.. math::
    \Delta t_\mathrm{a}=0.1\frac{\vert a\vert M_\mathrm{BH}}{s \dot{M}_\mathrm{BH,0}},

where :math:`s` is the spinup/spindown function. The numerical factor :math:`0.1` quantifies how finely we want to evolve spin; it ensures that the value of spin changes no more than :math:`10` per cent (relative to the current value) over the next time-step.

We also introduce a time-step related to the redirection of the spin vector. Since the spin vector may be redirected very quickly relative to its magnitude (due to LT torques), this criterion is separate to the one mentioned above. This time-step is given

.. math::
    \Delta t_\mathrm{a}=0.1\frac{M_\mathrm{warp}J_\mathrm{BH}}{\dot{M}_\mathrm{BH,0}J_\mathrm{warp}\sin\theta},

where :math:`\theta` is the angle between the current BH spin vector and the angular momentum of gas in the accretion disc on large scales. The numerical prefactor is again present to ensure a fine enough evolution of the spin vector direction. In particular, in the case that the spin vector and the gas angular momentum are perpendicular (:math:`\sin\theta=1`), this criterion will lead to a change of no more than :math:`\approx5\degree` in the spin vector direction per time-step.

