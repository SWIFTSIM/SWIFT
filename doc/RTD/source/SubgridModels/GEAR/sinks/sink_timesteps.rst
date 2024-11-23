.. Sink particles in GEAR model
   Darwin Roduit, 24 November 2024

Sink timesteps
~~~~~~~~~~~~~~

Sink particles interact with the surrounding gas through accretion. To accurately follow the local gas dynamics and assess the stability of gas hydrodynamics, we must impose timestep constraints on the sink particles.

First, sink particles obey a Courant-Friedrichs-Lewy (CFL)-like timestep constraint:

.. math::
   \Delta t_s \leq C_\text{CFL} \frac{r_{\text{cut, min}}}{\sqrt{c_{s,s}^2 + \| \Delta \mathbf{v}_s \|^2}} \; ,

where :math:`r_\text{cut, min} = \min(r_{\text{cut}, s}, \gamma_k \min_j(h_j))` is the minimal cut-off radius between the sink :math:`s` and the gas neihbours :math:`j`, :math:`c_{s, s}` and :math:`\Delta \mathbf{v}_s` are the gas soundspeed and the relative gas-sink velocity at the sink location. The latter two are reconstructed at the sink location with SPH interpolation.

With this condition, we ensure the sink knows the *local dynamics*. Since sink particles accrete gas, they must anticipate the surrounding gas infall. We achieve this with a gas free-fall time criterion, following Grudic 2017:

.. math::
   \Delta t_s \leq \sqrt{ \frac{3 \pi}{32 G \rho_s} } \quad \text{with} \quad \rho_s = \frac{3 m_s}{4 \pi {r_{\text{cut, min}}}} \; ,

with :math:`m_s` the sink mass.

These constraints ensure the gas accretion is smooth by removing few particles per timesteps. This is important since, in SPH, gas particles serve as interpolation points. By preventing the removal of huge amount of particles, we ensure the stability of the hydrodynamics computations.

To accurately follow *sink mergers*, we implemented Grudic 2017 two body timestep constraints between all sink particles:

.. math::
   \Delta t_s \leq \Delta t_\text{2-body} = \frac{t_\text{c, min} + t_\text{dyn, min}}{ t_\text{c, min} t_\text{dyn, min}} \; ,

with

.. math::
  \quad t_\text{c, min} = \min_{s \neq s'} \frac{\sqrt{ r_{ss'}^2 + \epsilon_{\text{sink}}^2} }{v_{ss'}} \quad \text{and} \quad t_\text{dyn, min} = \min_{s \neq s'} \sqrt{ \frac{(r_{ss'}^2 + \epsilon_\text{sink}^2)^{3/2}}{ G (m_s + m_{s'})}    } \; ,

where :math:`r_{ss'}` is the sinks relative separation, :math:`v_{ss'}` the relative velocity, :math:`m_{s}` and :math:`m_{s'}` their masses and :math:`\epsilon_{sink}` is the sink (fixed) gravitational softening.

We also implemented maximal timesteps sizes depending on the sink age; :math:`\Delta t_\text{max,s}^\text{age}`. A sink can be young, old or dead. In the first two cases, the sink's timestep is :math:`\min(\Delta t_\text{max,s}^\text{age}, \Delta t_s)`. In the last case, we impose :math:`\Delta t_\text{2-body}` only if a dead sink is involved in a two-boy encounter with an alive sink. Otherwise, the sink has no timesteps constraint (apart from gravity).

Notice that sink particles also satisfy a gravity timestep constraint, as all gravitational particles in Swift. We add a last timestep contraint, explained in the following section.


Star formation constraint
=========================

A concern that we did not adressed know concerns the star spawning algorihtm. Althought the stars' masses are sampled from the IMF, the stars metallicities are not. If a sink accrete lots of mass, then it can create lots of star particles during the same timestep. However, these stars will all have the same metallicity, which does not represent the actual metals evolution during the accretion. In such situation, the galaxies' properties are affected and do not represent the underlying physics. Another problem is that we can spawn many stars at once and the code may complain. Such situation is not ideal. Although, the constraints in the previous section will help, they are not sufficient. Hence, we introduce a new timestep constraint.
