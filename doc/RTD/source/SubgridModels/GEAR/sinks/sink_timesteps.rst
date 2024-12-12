.. Sink particles in GEAR model
   Darwin Roduit, 24 November 2024

.. _sink_GEAR_timesteps:

Sink timesteps
~~~~~~~~~~~~~~

Sink particles interact with the surrounding gas through accretion. To accurately follow the local gas dynamics and assess the stability of gas hydrodynamics, we must impose timestep constraints on the sink particles.

First, we want to ensure sinks know the *local dynamics* and thus they obey a Courant-Friedrichs-Lewy (CFL)-like timestep constraint:

.. math::
   \Delta t_s \leq \Delta t_\text{CFL} =  C_\text{CFL} \frac{r_{\text{cut, min}}}{\sqrt{c_{s,s}^2 + \| \Delta \mathbf{v}_s \|^2}} \; ,

where :math:`r_\text{cut, min} = \min(r_{\text{cut}, s}, \gamma_k \min_j(h_j))` is the minimal cut-off radius between the sink :math:`s` and the gas neighbours:math:`j`, :math:`c_{s, s}` and :math:`\Delta \mathbf{v}_s` are the gas sound-speed and the relative gas-sink velocity at the sink location. The latter two are reconstructed at the sink location with SPH interpolation. The value of :math:`C_\text{CFL}` is given in the YAML parameter file with ``GEARSink:CFL_condition``.

Since sink particles accrete gas, they must anticipate the surrounding gas infall. We achieve this with a gas free-fall time criterion similar to `Grudic et al. (2021) <https://academic.oup.com/mnras/article/506/2/2199/6276745>`_:

.. math::
   \Delta t_s \leq \Delta t_\text{ff} = \sqrt{ \frac{3 \pi}{32 G \rho_s} } \quad \text{with} \quad \rho_s = \frac{3 m_s}{4 \pi {r_{\text{cut, min}}}} \; ,

with :math:`m_s` the sink mass.

These constraints ensure smooth gas accretion by removing a few particles per timestep. This is important since, in SPH, gas particles serve as interpolation points. By preventing the removal of a large number of particles, we ensure the stability of the hydrodynamic computations.

Sink 2-body encounters
++++++++++++++++++++++

To accurately follow *sink mergers*, we implemented `Grudic et al. (2021) <https://academic.oup.com/mnras/article/506/2/2199/6276745>`_ two body timestep constraints between all sink particles:

.. math::
   \Delta t_s \leq \Delta t_\text{2-body} = \frac{ t_\text{c, min} t_\text{dyn, min}}{t_\text{c, min} + t_\text{dyn, min}} \; ,

with

.. math::
  \quad t_\text{c, min} = \min_{s \neq s'} \frac{ |\varphi^{-1}(r_{ss'}, \, H_\text{sink})| }{v_{ss'}} \quad \text{and} \quad t_\text{dyn, min} = \min_{s \neq s'} \sqrt{ \frac{ |\varphi'(r_{ss'}, \, H_\text{sink})|^{-1}} { G (m_s + m_{s'})}    } \; ,

where :math:`r_{ss'}` is the sinks relative separation, :math:`v_{ss'}` the relative velocity, :math:`m_{s}` and :math:`m_{s'}` their masses and :math:`H_\text{sink}` is the sink (fixed) gravitational softening. The function :math:`\varphi(r, H)` is the potential corresponding to the Wendland C2 kernel density field (see `Schaller et al. (2024) <https://doi.org/10.1093/mnras/stae922>`_ section 4.1) and  :math:`\varphi'(r, H) \equiv \frac{\mathrm{d} \varphi(r, H)}{\mathrm{d} r}` its derivative.

Timesteps per sink's age categories
+++++++++++++++++++++++++++++++++++

We also implemented maximal timesteps sizes depending on the sink age; :math:`\Delta t_\text{max,s}^\text{age}`. A sink can be young, old or dead. In the first two cases, the sink's timestep is :math:`\min(\Delta t_\text{max,s}^\text{age}, \Delta t_s)`. In the last case, we impose :math:`\Delta t_\text{2-body}` only if a dead sink is involved in a two-boy encounter with an alive sink. Otherwise, the sink has no timestep constraint (apart from gravity). The parameters controlling the transition between the young and old is ``GEARSink:timestep_age_threshold``, and the one between old and dead is ``GEARSink:timestep_age_threshold_unlimited_Myr``. The maximal timesteps are given by  ``GEARSink:max_timestep_young_Myr`` and  ``GEARSink:max_timestep_old_Myr``.

Notice that sink particles also satisfy a gravity timestep constraint, as do all gravitational particles in Swift.

Star formation constraint
+++++++++++++++++++++++++

Although the stars' masses are sampled from the IMF, the stars' metallicities are not. If a sink accretes mass, it can create many star particles simultaneously. However, these stars will all have the same metallicity, which does not represent the actual metals' evolution during the accretion. In such a situation, the galaxies' properties are affected and do not represent the underlying physics.

Another problem is that we can spawn many stars simultaneously, and the code may complain. Such a situation could be better. Although the constraints in the previous section will help, more is needed. Our solution is to introduce a new accretion criterion using the IMF properties. However, since our politics is that accretion should be feedback-regulated and not based on an arbitrary accretion rate, we reduce the sink time step to avoid limiting the star formation rate to an arbitrary value.

The new accretion criterion is the following. The swallowed gas and sink mass does not exceed ``n_IMF`` times the IMF mass (see the IMF sampling section), but make sure to swallow at least one particle: :math:`M_\text{swallowed} \leq n_\text{IMF} M_\text{IMF} \text{ or } M_\text{swallowed} = 0`.

Since we artificially restrict mass accretion, we keep track of the mass :math:`M_\text{eligible}` that would be swallowed without this criterion. Then, we compute the error :math:`\Delta M` between the restricted and unrestricted swallow. The absolute error is :math:`\Delta M = M_\text{swallowed} - M_\text{eligible}` and the relative error is :math:`| \Delta M | / M_\text{eligible}`.

When :math:`\Delta M < 0` (i.e. :math:`M_\text{swallowed} \neq M_\text{eligible}`), we know the accretion was restricted and we can apply another time-step contraint. To compute a timestep, we convert :math:`\Delta M` to accretion rate by dividing by :math:`\Delta t_\text{s, estimated} = \min(\Delta t_\text{CFL}, \, \Delta t_\text{ff}, \Delta  t_\text{2-body})`. Hence, we have the constraint:

.. math::
   \Delta t_s \leq \Delta t_\text{SF} = \eta \cfrac{M_\text{eligible} \Delta t_\text{s, estimated}}{\Delta M} \text{ if } \Delta M < 0 \; ,

where :math:`\eta` is a tolerance parameter. This parameter corresponds to ``GEARSink:tolerance_SF_timestep`` in the code.

