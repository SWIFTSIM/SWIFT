.. AGN spin and jet model
   Filip Husko, 1 April 2022

.. AGN_spin_jet:

Jet launching algorithm
-----------------------

In order to launch jets, we introduce a jet reservoir that functions identically to the thermal reservoir used in EAGLE and COLIBRE. When the jet reservoir exceeds the value :math:`N_\mathrm{j}\overline{m}_\mathrm{ngb}v_\mathrm{j}^2/2`, kicks are given out to :math:`N_\mathrm{j}` particles. Here, :math:`N_\mathrm{j}` is the target number of particles to kick per each kicking event (typically equal to 2, but should always be a multiple of 2 to ensure approximately symmetrical jets). :math:`\overline{m}_\mathrm{ngb}v_\mathrm{j}^2/2` is the average kinetic energy of one kicking event, with :math:`\overline{m}_\mathrm{ngb}` the mean neighbour mass of gas particles in the BH smoothing kernel and :math:`v_\mathrm{j}` the target kicking velocity.

These kicks are handed out to particles in a symmetric way with respect to the spin vector of the BH. :math:`N_\mathrm{j}/2` particles are kicked from the 'upper' hemisphere relative to the spin vector, and the other half from the lower hemisphere. The particles to be kicked can be any in the smoothing kernel. We include four different choices: the particles kicked are: 1) the closest to the BH, 2) the farthest from the BH, 3) the ones of minimal density and 4) the ones closest to the spin axis, in terms of angular distance. Note that these sortings are done for each hemisphere seperately. 

The particles chosen are always given velocities based on the same algorithm, regardless of their positions in the kernel. We perform the actual kicks in the following way. Velocity kicks are chosen at random from a cone around the current spin vector with a (half-)opening angle of :math:`\theta_\mathrm{j}`. In particular, we first choose the kick vector around the z-axis as :math:`\vec{v}_\mathrm{kick}=(\sin\theta\cos\phi,\hspace{0.3mm}\sin\theta\sin\phi,\hspace{0.3mm}\cos \theta)`. Here, :math:`\cos\theta` is chosen uniformly from the interval :math:`[\cos\theta_\mathrm{j},1]`, and :math:`\sin\theta=\sqrt{1-\cos\theta^2}`. :math:`\phi` is chosen uniformly from :math:`[0,2\pi]`. This random vector, now representing a random kick within a cone around the z-axis, is rotated into the frame of the spin vector so that the cone is pointing in the right direction. For particles being kicked from the 'negative' side of the BH hemisphere, the final kick vector is simply multiplied by :math:`-1`.

We then add the kick vector to the particle's current velocity. We do this in a way that conserves energy, so that the magnitude of the final velocity is computed from

.. math::
    \frac{1}{2}m_\mathrm{gas}\vec{v}_\mathrm{fin}^2=\frac{1}{2}m_\mathrm{gas}\vec{v}_\mathrm{0}^2 + \frac{1}{2}m_\mathrm{gas}\vec{v}_\mathrm{kick}^2,
    
while its direction is computed from conservation of momentum:

.. math::
    m_\mathrm{gas}\vec{v}_\mathrm{fin}=m_\mathrm{gas}\vec{v}_\mathrm{0} + m_\mathrm{gas}\vec{v}_\mathrm{kick}.

Black hole time steps
---------------------

Black holes will generally have time steps based on their gravitational interactions, but also based on their current accretion rate and the expected time interval until which the thermal feedback reservoir (representing radiative feedback) will have grown enough to heat 1 particle. We introduce a similar time step, but based on when the jet reservoir grows enough to kick :math:`N_\mathrm{j}` particles. We then take the minimum of those two. 

On top of that, we add a time step that makes sure the BH spin doesn't get evolved too much over one time step. Its magnitude is actually not a problem, since the growth of the spin magnitude is always tied with the growth of mass. However, the direction is more problematic (especially in the thin disk case, where alignment can occur very quickly, with little mass or spin growth, due to large warp radii). For this reason, we make sure that the amount of warp angular momentum interacting with the BH over the next time-step, :math:`\Delta_J=(J_\mathrm{warp}/M_\mathrm{warp})\Delta M`, is a small fraction of the current BH angular momentum :math:`J_\mathrm{BH}` (e.g. :math:`0.01`).

