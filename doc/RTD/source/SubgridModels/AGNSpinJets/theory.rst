.. AGN spin and jet model
   Filip Husko, 1 April 2022

.. AGN_spin_jet:


Model summary
-------------

Here we provide a comprehensive summary of the model. In order to more easily visualize the model, it is recommended to run the python script in order to generate plots that show various aspects of the model.

Any model for realistic AGN jets must include black hole spin since jet powers depend steeply on spin, and because including spin provides a well-defined direction for the jets to be launched in. The spin (angular momentum) of BHs is best represented through the dimensionlesss spin parameter :math:`a=J_\mathrm{BH}c/M_\mathrm{BH}^2 G`, where :math:`J_\mathrm{BH}` is its angular momentum. For theoretical reasons, its magnitude cannot grow above 1. It can be positive, representing prograde accretion, or negative, representing retrograde accretion.

Jet powers, in addition to depending on spin, also depend on which accretion state the black hole is in. We refer to these states by the shape of the accretion disk that surrounds the BH. We include three accretion states: the thick (or advection-dominated accretion flow; ADAF), thin (standard) and slim disk. Our main reference points for these disks are the following papers: `Shakura & Sunyaev (1973) <https://ui.adsabs.harvard.edu/abs/1973A%26A....24..337S/abstract>`_, `Narayan & Yi (1994) <https://ui.adsabs.harvard.edu/abs/1994ApJ...428L..13N/abstract>`_ and `Wang & Zhou. (1999) <https://ui.adsabs.harvard.edu/abs/1999ApJ...516..420W/abstract>`_, respectively.

The thick disk appears at low accretion rates, has very strong jets and is inefficient at spinning up the black hole. The thin disk, appearing at intermediate accretion rates, typically has weak jets, strong radiation and efficiently spins up the black hole. The slim disk, corresponding to super-Eddington accretion, has features of both: in terms of geometry, orbits and angular momentum, it is similar to the thick disk. It is optically thin, leading to strong radiation. However, it also has strong jets. We assume that each subgrid accretion disk launches jets and radiates at the same time, regardless of the type it is. However, we use expressions for the jet and radiative efficiencies that depend on the type of the disk, and which are physically motivated.

Transitions from one accretion state to another
-----------------------------------------------

.. figure:: modes.png
    :width: 600px
    :align: center
    :figclass: align-center
    :alt: Accretion regimes

    The type of accretion disk surrounding the BHs depending on their accretion rates and spins. The transition between the thick and thin disk is calculated assuming the viscosity parameter :math:`\alpha=0.2`, while the transition from thin to slim disk is assumed to occur when the latter is :math:`F=0.5` times as radiatively efficienct as the former.

The state of the subgrid accretion disk depends mostly on the Eddington fraction, i.e. the (dimensionless) accretion rate of the BH in units of the Eddington accretion rate, which we denote as :math:`\dot{m}`. We assume that the subgrid accretion disk is thick for :math:`\dot{m}<0.03`, based on observations (`Russell et al. 2013 <https://ui.adsabs.harvard.edu/abs/2013MNRAS.432..530R/abstract>`_). This also allows us to constrain the value of one of the main parameters in any accretion model: the viscosity parameter :math:`\alpha` (which is related to the kinematic viscosity :math:`\nu`and sound speed :math:`c_\mathrm{s}` through :math:`\nu=\alpha c_\mathrm{s}H`, with :math:`c_\mathrm{s}` the sound speed and :math:`H` the disk half-thickness). Numerical calculations suggest that thick disks are present for :math:`\dot{m}<0.4\alpha^2` (`Yuan & Narayan 2014 <https://ui.adsabs.harvard.edu/abs/2014ARA%26A..52..529Y/abstract>`_), and this agrees with observations if :math:`\alpha=0.25-0.3`. These values agree very well with more direct observational estimates, which suggest :math:`\alpha=0.2-0.3` (`Martin et al. 2019 <https://ui.adsabs.harvard.edu/abs/2014ARA%26A..52..529Y/abstract>`_).

The transition from the thin to the slim disk should occur around :math:`\dot{m}\approx 1`. However, the exact physics of this transition is not well understood. There is likely some spin dependence of the critical accretion rate, due to different radiative physics depending on spin. One of the main properties of slim disks is that they are less radiatively efficient than thin disks (`Sadowski et al. 2014 <https://ui.adsabs.harvard.edu/abs/2014MNRAS.439..503S/abstract>`_). We thus assume that the transition occurs when the radiative efficiency of a slim disk, :math:`\epsilon_\mathrm{r,SD}`, falls below some fraction of the radiative efficiency of a thin disk, :math:`\epsilon_\mathrm{r,TD}`. We quantify this as :math:`\epsilon_\mathrm{r,SD}<F\epsilon_\mathrm{r,SD}`, with :math:`F\approx 0.5` a free parameter. We give the expressions for both of the efficiencies below.

Jet efficiencies
----------------

The jet efficiency is related to the jet power through :math:`\epsilon_\mathrm{j}=P_\mathrm{j}/\dot{M}_\mathrm{BH,0}c^2`, where :math:`\dot{M}_\mathrm{BH,0}` is the accretion rate measured in the simulation, e.g. the Bondi rate). We use the formula for the jet efficiency based on general-relativistic, magneto-hydrodynamical (GRMHD) simulations by `Tchekhovskoy et al. (2010) <https://ui.adsabs.harvard.edu/abs/2010ApJ...711...50T/abstract>`_:

.. math::
    \epsilon_\mathrm{j}=\frac{\kappa}{4\pi}\bigg(\frac{H/R}{0.3}\bigg)^\eta \phi_\mathrm{BH}^2\Omega_\mathrm{BH}^2\big(1+1.38\Omega_\mathrm{BH}^2-9.2\Omega_\mathrm{BH}^4\big),

where :math:`\kappa\approx0.05` is a numerical factor which depends on the initial geometry of the magnetic field, :math:`\phi_\mathrm{BH}` is the dimensionless magnetic flux threading the horizon (see original paper for precise definition), and :math:`\Omega_\mathrm{BH}=a/2r_\mathrm{H}` is the (dimensionless) angular velocity of the black hole event horizon. Here, :math:`r_\mathrm{H}=1+\sqrt{1-a^2}` is the radius of the horizon in units of the gravitational radius :math:`R_\mathrm{G}=M_\mathrm{BH}G/c^2`. The formula above, for the jet efficiency, agrees very well with the results from higher-resolution simulations performed by `Narayan et al. (2021) <https://ui.adsabs.harvard.edu/abs/2010ApJ...711...50T/abstract>`_, who provide the following fit for the magnetic flux as a function of spin:

.. math::
    \phi_\mathrm{BH}(a)=-20.2a^3-14.9a^2+34a+52.6.
    
The `Tchekhovskoy et al. (2010) <https://ui.adsabs.harvard.edu/abs/2010ApJ...711...50T/abstract>`_ jet efficiency depends very steeply on spin (:math:`\epsilon_\mathrm{j}\propto a^2` for small spin and :math:`\epsilon_\mathrm{j}\propto a^6` near :math:`a=1`). It can reach values above 100 per cent for large spins, and is also different (weaker) for negative spins.

The dependence of the jet efficiency on the type of accretion disk comes through the factor that depends on the aspect ratio :math:`H/R`, since accretion disks differ in this quantity. Theoretical, self-similar models of thick disks suggest :math:`H/R=0.5` (`Narayan & Yi 1995b <https://ui.adsabs.harvard.edu/abs/1995ApJ...452..710N/abstract>`_), but we instead take :math:`H/R=0.3`, more in line with simulations. For slim disks, which have received less attention in simulations, we assume the value :math:`H/R=1/(2\sqrt{5})\approx 0.2` (based on the self-similar model by `Wang & Zhou. 1999 <https://ui.adsabs.harvard.edu/abs/1999ApJ...516..420W/abstract>`_).

Thin disks are, not surprisingly, much thinner. The value of :math:`H/R` in this regime is not a constant, but rather depends on the BH mass and accretion rate, slightly on radius and also on the viscosity parameter :math:`\alpha`. Thin disks have three different regions in the `Shakura & Sunyaev (1973) <https://ui.adsabs.harvard.edu/abs/1973A%26A....24..337S/abstract>`_ model. For simplicity, we model the whole disk as being represented with only one region. In region a), the innermost one, radiation dominates over gas pressure. It is typically very small or doesn't exist at all, so we disregard it as a possibility. In regions b) and c), gas pressure dominates over radiation pressure. In b), electrons dominate in the opacity, while in c), free-free absorption dominates. We leave both regions as a possibility, and leave the choice as a free parameter in the model (not likely to lead to large differences in galaxy/BH evolution). The expressions for the aspect ratio in these regions are

.. math::
    \bigg(\frac{H}{R}\bigg)_\mathrm{TD,b} = 1.25\times10^{-3} \alpha^{-1/10}\dot{m}^{1/5}\bigg(\frac{M_\mathrm{BH}}{10^8\hspace{0.5mm}\mathrm{M}_\odot}\bigg)^{-1/10}\bigg(\frac{R}{2R_\mathrm{G}}\bigg)^{1/20}

in region b) and

.. math::
    \bigg(\frac{H}{R}\bigg)_\mathrm{TD,c} = 1.15\times10^{-3} \alpha^{-1/10}\dot{m}^{3/20}\bigg(\frac{M_\mathrm{BH}}{10^8\hspace{0.5mm}\mathrm{M}_\odot}\bigg)^{-1/10}\bigg(\frac{R}{2R_\mathrm{G}}\bigg)^{1/8}

in region c). 

The jet efficiency also depends on the slope :math:`\eta`. Classical jet theory (`Meier 2001 <https://ui.adsabs.harvard.edu/abs/2001Sci...291...84M/abstract>`_) suggests that jet powers depend on the aspect ratio linearly, so :math:`\eta=1`. This is also in line with some simulations finding a reduction in jet efficiencies with the aspect ratio (e.g. `Tchekhovskoy et al. 2014 <https://ui.adsabs.harvard.edu/abs/2014MNRAS.437.2744T/abstract>`_). In this scenario, jets launched from thin disks are of order :math:`\approx100` times less powerful than those launched from thick disks. On the other hand, some simulations of thin disks have found jet efficiencies similar to thick disk ones (e.g. `Liska et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..550L/abstract>`_), which is supported by observations of blazars (`Ghisellini et al. 2014 <https://ui.adsabs.harvard.edu/abs/2014Natur.515..376G/abstract>`_). In this picture, thin jets are approximately as efficient as thick disk ones, which can in our case be implemented as :math:`\eta=0`. The reality is likely to be somwhere in between. Note that the choice of :math:`\eta` likely has a strong impact on the evolution of galaxies and BHs; our default choice is the classical picture in which :math:`\eta=1`.

.. figure:: efficiencies.png
    :width: 1200px
    :align: center
    :figclass: align-center
    :alt: Efficiencies

    Feedback efficiencies (jet - blue, radiation - red) for all three accretion disk types. Shaded regions represent likely ranges of efficiencies (where the efficiencies depend on mass and/or accretion rate). The thin disk jet efficiencies were computed assuming the slope of the efficiency vs. aspect ratio relation is :math:`\eta=1`, and the aspect ratios were computed for region b) of the Shakura & Sunyaev solution. Radiative efficiencies in the thick disk were computed assuming the electron heating parameter :math:`\delta=0.2`.

Radiative efficiencies
----------------------

In the EAGLE and COLIBRE models, all subgrid accretion disks are effectively thin, and the BH is always assumed to be in this regime. In our model, the radiative efficiency (defined in an analagous way to the jet efficiency, but using the luminosity) is no longer fixed at a value of order :math:`10` per cent. Instead, we use spin-dependant formulas that vary with the type of disk. In the thin disk, the radiative efficiency :math:`\epsilon_\mathrm{r,TD}` is related to the binding energy at the innermost stable circular orbit (ISCO) and is given by

.. math::
    \epsilon_\mathrm{r,TD}(a) = 1-e_\mathrm{ISCO}(a)=1-\sqrt{1-\frac{2}{3r_\mathrm{ISCO}(a)}}.
    
Here, :math:`r_\mathrm{ISCO}` is the radius of the ISCO in gravitational radii (see e.g. appendix B of `Fiacconi et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.3807F/abstract>`_ for an expression giving the spin dependence). The radiative efficiency of the thin disk grows slowly from its minimum value of :math:`\approx4` per cent for :math:`a=-1` to :math:`\approx5.5` per cent for :math:`a=0`. For positive spins it grows more steeply; it is :math:`10` per cent by :math:`a=0.65`. Beyond that the dependence steepens even further, with values of :math:`20`, :math:`30` and :math:`40` per cent reached at :math:`a=0.95`, :math:`a=0.997` and :math:`a=1`, respectively.

In the thick disk regime, radiative efficiencies are lower by a factor :math:`\approx100` than jet efficiencies. The formulas we use are based on results by `Mahadevan (1997) <https://ui.adsabs.harvard.edu/abs/1997ApJ...477..585M/abstract>`_, who studied cooling processes of electrons (which dominate in the radiation) in the context of the original thick disc solution. They found two different regimes: for :math:`\dot{m}<\dot{m}_\mathrm{crit,visc}`, viscous heating dominates the heating of electrons, whereas for :math:`\dot{m}_\mathrm{crit,visc}<\dot{m}<\dot{m}_\mathrm{crit,ADAF}`, it is dominated by ion-electron heating. Here, :math:`\dot{m}_\mathrm{crit,visc}` is the transitional value between the two thick disc (ADAF) regimes, and :math:`\dot{m}_\mathrm{crit,ADAF}=0.4\alpha^2` is the transitional accretion rate which separates thin and thick discs. The radiative efficiency in the viscous heating regime is given by

.. math::
    \epsilon_\mathrm{r,ADAF}=0.0002\epsilon_\mathrm{r,TD}\bigg(\frac{\delta_\mathrm{ADAF}}{0.0005}\bigg)\bigg(\frac{1-\beta}{0.5}\bigg)\bigg(\frac{6}{r_\mathrm{ISCO}}\bigg),

while in the ion-heating regime it is given by

.. math::
    \epsilon_\mathrm{r,ADAF}=0.2\epsilon_\mathrm{r,TD}\bigg(\frac{\dot{m}}{\alpha^2}\bigg)\bigg(\frac{\beta}{0.5}\bigg)\bigg(\frac{6}{r_\mathrm{ISCO}}\bigg).
    
Here, :math:`\beta` is the ratio of gas pressure and total pressure (which includes the magnetic pressure). `Yuan & Narayan (2014) <https://ui.adsabs.harvard.edu/abs/2014ARA%26A..52..529Y/abstract>`_ define a somewhat different parameter, :math:`\beta_\mathrm{ADAF}`, as the ratio of gas pressure and magnetic pressure. The two parameters are related by :math:`\beta=\beta_\mathrm{ADAF}/(1+\beta_\mathrm{ADAF})`. :math:`\beta_\mathrm{ADAF}` is not an independent parameter; many simulations have found that :math:`\alpha\beta_\mathrm{ADAF}\approx0.5` (e.g. `Begelman et al. 2021 <https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.2040B/abstract>`_, see also `Yuan & Narayan 2014 <https://ui.adsabs.harvard.edu/abs/2014ARA%26A..52..529Y/abstract>`_ for a review), which we adopt. :math:`\delta_\mathrm{ADAF}` represents the fraction of viscous energy transferred to the electrons, and is constrained in theoretical studies between 0.1 and 0.5 (`Yuan & Narayan 2014 <https://ui.adsabs.harvard.edu/abs/2014ARA%26A..52..529Y/abstract>`_, `Sharma et al. 2007 <https://ui.adsabs.harvard.edu/abs/2007ApJ...667..714S/abstract>`_). Observations imply a value close to 0.2 (`Yuan et al. 2003 <https://ui.adsabs.harvard.edu/abs/2003ApJ...598..301Y/abstract>`_, `Liu & Wu 2013 <https://ui.adsabs.harvard.edu/abs/2013ApJ...764...17L/abstract>`_). The critical accretion rate between the two thick disc regimes can be found by ensuring that both formulas presented above yield the same radiative efficiency (at that accretion rate). This gives an accretion rate equal to

.. math::
    \dot{m}_\mathrm{crit,visc}=0.0002\bigg(\frac{\delta_\mathrm{ADAF}}{0.0005}\bigg)\bigg(\frac{1-\beta}{\beta}\bigg)\alpha^2.
    
For slim disks we take the radiative efficiency based on GRMHD simulations of super-Eddington accretion (for various BH spins) performed by `Sadowski et al. (2014) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.439..503S/abstract>`_. `Madau et al. (2014) <https://ui.adsabs.harvard.edu/abs/2014ApJ...784L..38M/abstract>`_ found the following fitting function which represents the `Sadowski et al. (2014) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.439..503S/abstract>`_ results:

.. math::
    \epsilon_\mathrm{r,SD}=\frac{0.1}{\dot{m}}A(a)\bigg( \frac{0.985}{1.6/\dot{m}+B(a)}+\frac{0.015}{1.6/\dot{m}+C(a)}\bigg),
    
where the three spin-dependant functions are given by :math:`A(a)=(0.9663-0.9292a)^{-0.5639}`, :math:`B(a)=(4.627-4.445a)^{-0.5524}` and :math:`C(a)=(827.3-718.1a)^{-0.7060}`. The radiative efficiency of slim disks, based on this formula, matches the thin disk radiative efficiency (given at the beginning of the section) at low accretion rates. At high accretion rates (:math:`\dot{m}\gtrapprox1`, but depending on spin), the radiative efficiency drops. These two formulas are used to decide when a disk transitions from thin to slim.

Evolution of the black hole spin magnitude
------------------------------------------

.. figure:: spec_ang_mom.png
    :width: 600px
    :align: center
    :figclass: align-center
    :alt: Angular momenta

    Dimensionless pecific angular momentum of the thin disk at the innermost stable circular orbit (ISCO, solid red line), compared with the specific angular momentum at the inner radius (the event horizon) for advection-dominated flows (the thick and slim disk) for a few values of the viscosity parameter :math:`\alpha`. The dashed red line shows that the latter can be approximated as :math:`45` per cent of the former.

The BH spin (or angular momentum) is, naturally, a vector. However, due to Lense-thirring torques (we discuss these in more detail below), the accretion disk is always either aligned or counteraligned with the rotational axis of the black hole. This means that almost all relevant quantities, such as the efficiencies discussed above, can be expressed as depending only on the magnitude of spin, but also allowing for a negative sign to account for counteraligned disks (retrograde accretion). This is also true for the evolution of the magnitude of spin.

In the absence of jet spindown, the evolution of angular momentum is given simply by :math:`\dot{J}_\mathrm{BH}=L_\mathrm{in}\dot{M}_\mathrm{BH}`, where :math:`L_\mathrm{in}` is the specific angular momentum at the inner radius of the accretion disk. This can be transformed into an equation for spin evolution, yielding

.. math::
    \frac{\mathrm{d}a}{\mathrm{d}\ln M_\mathrm{BH,0}}=\ell_\mathrm{in}-2a e_\mathrm{in},
    
where :math:`\ell_\mathrm{in}` is the specific angular momentum in units where :math:`G` and :math:`c` are equal to unity, and :math:`\mathrm{d}\ln M_\mathrm{BH,0}=\mathrm{d}M_\mathrm{BH,0}/M_\mathrm{BH}` is the logarithmic change in mass, not including losses due to radiation (`Fanidakis et al. 2011 <https://ui.adsabs.harvard.edu/abs/2011MNRAS.410...53F/abstract>`_). The specific binding energy can be related to the radiative efficiency through :math:`e_\mathrm{in}=1-\epsilon_\mathrm{r}` for all three accretion states (for the thick disc, the radiative efficiency is negligible for this application). 

For the thin disc, the inner radius :math:`R_\mathrm{in}` can be taken to be the radius of the ISCO, since orbits should quickly degrade within it. We thus take :math:`\ell_\mathrm{in}` as the specific angular momentum at the ISCO for the thin disc (the expression for which is given in e.g. Appendix B of `Fiacconi et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.3807F/abstract>`_). For the thin disk, the spinup function (the equation shown above) is always positive, meaning that the BH will always be spun up to positive values. This means that the BH will be spun down if spin is negative, or spun up to an equilibrium value of :math:`a_\mathrm{eq}=1` if spin is positive. For advection-dominated disks (the thick and slim disk), we assume that :math:`\ell_\mathrm{in}` is :math:`45` per cent of the ISCO value, based on numerical GR calculations by `Popham & Gammie (1998) <https://ui.adsabs.harvard.edu/abs/1998ApJ...504..419P/abstract>`_. We base this assumption on fits of the `Popham & Gammie (1998) <https://ui.adsabs.harvard.edu/abs/1998ApJ...504..419P/abstract>`_ results done by `Benson & Babul (2009) <https://ui.adsabs.harvard.edu/abs/2009MNRAS.397.1302B/abstract>`_. We independently compared these fits to the ISCO values and found :math:`\ell_\mathrm{in}\approx0.45\ell_\mathrm{ISCO}` with no more than :math:`10` per cent error for all values of spin and relevant values of :math:`\alpha=0.1-0.3`.

For the thick and slim disk, these lower specific angular momenta lead to a non-zero equilibrium spin value :math:`a_\mathrm{eq}<1`. If :math:`a>a_\mathrm{eq}`, the BH will be spun down due to frame-dragging and viscosity; the frame-dragging rotationally accelerates any accreting gas (on account of the BH angular momentum), while viscosity carries away some of that angular momentum. Including jets into the model leads to further spindown. The jet spindown term (to be added to the spinup equation above) can be derived as 

.. math::
    \bigg(\frac{\mathrm{d}a}{\mathrm{d}\ln M_\mathrm{BH,0}}\bigg)_\mathrm{j}=-\epsilon_\mathrm{j}(a)\frac{\sqrt{1-a^2}}{a}\bigg[\Big(\sqrt{1-a^2}+1 \Big)^2+a^2 \bigg]
    
(see `Benson & Babul 2009 <https://ui.adsabs.harvard.edu/abs/2009MNRAS.397.1302B/abstract>`_ for a derivation, which we have independently verified). Including jet spindown leads to even lower equilibrium spin values; e.g. for the thick disk this is only :math:`a_\mathrm{eq}\approx0.25`.

.. figure:: spinup.png
    :width: 1200px
    :align: center
    :figclass: align-center
    :alt: Spinup/spindown function

    Spinup/spindown function (the rate of black hole spin evolution) as a function of spin for all three accretion disk types. Black lines show evolution with only accretion included, while blue lines show the total including jet spindown. These plots show that the thin disk is always spun up to to :math:`a_\mathrm{eq}=1`, even with jets (due to low jet efficiencies). The advection-dominated disks (thick and slim disk) are spun up to positive equilibrium values :math:`a_\mathrm{eq}<1`, or spun down to such an equilibrium value if :math:`a>a_\mathrm{eq}`. This is due to extraction of rotational energy from the BH by frame dragging and transport of the angular momentum to large distances through viscous forces. Including jet spindown pushes these equilibrium spins to even smaller values.

Evolution of the black hole spin direction
------------------------------------------

In the previous section we claimed that the evolution of the magnitude of spin depends only on whether accretion is prograde or retrograde. The two remaining questions are: 1) what about its direction, and 2) how to decide whether accretion is prograde or retrograde. We will now address the first of these.

Lense-Thirring torques (`Lense & Thirring 1918 <https://ui.adsabs.harvard.edu/abs/1918PhyZ...19..156L/abstract>`_) arise from additional GR forces that operate near spinning BHs, related to the frame-dragging of spacetime. In isolation, they cause the precession of a parcel of gas as it orbits around the BH. For accretion disks, their effects depend on the type of disk (see `Nixon & King 2016 <https://ui.adsabs.harvard.edu/abs/2016LNP...905...45N/abstract>`_ for a review). Lense-Thirring torques do not have a component in the direction of the BH spin vector, which is why they do not play a role in the evolution of the magnitude of spin.

In all cases, Lense-Thirring torques are effective only within some radius :math:`R_\mathrm{warp}`, which marks the boundary between the outer disk and an inner region, within which the BH can 'communicate' through these torques with the disk. Within this radius, the disk is on average aligned or counteraligned with the BH, whereas outside it, it is aligned with some large-scale angular momentum direction (which we can measure in the simulation) - hence the name warp radius. Given some surface density, one can also define the warp mass :math:`M_\mathrm{warp}` and the warp angular momentum :math:`J_\mathrm{warp}` as the total mass and angular momentum within :math:`R_\mathrm{warp}`, respectively. We will discuss how all of these warp-related quantities are calculated in each of the accretion disks further below, but for now we focus on how these warped disks feature in our model.

In terms of the evolution of the spin direction, the main assumption of our model is as follows (see `King et al. 2005 <https://ui.adsabs.harvard.edu/abs/2005MNRAS.363...49K/abstract>`_ for the original argument, and additional discussions in e.g. `Fanidakis et al. 2011 <https://ui.adsabs.harvard.edu/abs/2011MNRAS.410...53F/abstract>`_, `Fiacconi et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.3807F/abstract>`_ and `Griffin et al. 2019a <https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..198G/abstract>`_). All matter that flows through an accretion disk is aligned or counteraligned with the BH spin vector in the accretion process. Due to conservation of angular momentum, the spin vector itself also has to adjust to keep the total angular momentum conserved. In the process of consuming one warp mass :math:`M_\mathrm{warp}`, the direction of the BH spin vector is aligned to match the direction of the total angular momentum of the system comprising the BH and the disk out to the warp radius. The direction of the BH spin vector can then be determined from :math:`\vec{J}_\mathrm{warp}=\vec{J}_\mathrm{BH}+J_\mathrm{warp}\hat{J}_\mathrm{d}`, where :math:`\vec{J}_\mathrm{BH}` is the old BH angular momentum vector, and :math:`\hat{J}_\mathrm{d}` is the direction of the large-scale accretion disk (which we assume matches the direction of the angular momentum of the gas in the BH smoothing kernel).

In practice, the BH will consume parcels of mass that differ from :math:`M_\mathrm{warp}`. We assume that any such parcel of mass :math:`\Delta M` (e.g. the mass to be consumed within a single time step) can be split up onto :math:`n=\Delta M / M_\mathrm{warp}` individual increments of accretion, so the total angular momentum of the system within that time step is :math:`\vec{J}_\mathrm{warp}=\vec{J}_\mathrm{BH}+n J_\mathrm{warp}\hat{J}_\mathrm{d}`, i.e. :math:`n` warp angular momenta are consumed, with an angular momentum of :math:`\Delta \vec{J}=n J_\mathrm{warp}\hat{J}_\mathrm{d}=(J_\mathrm{warp}/M_\mathrm{warp})\Delta M `. This can also be viewed as the BH consuming material with a specific angular momentum of :math:`L_\mathrm{warp}=J_\mathrm{warp}/M_\mathrm{warp}`. Note that this picture is only valid if the BH spin vector does not change much during this process (in both magnitude and direction), which can be ensured with wisely chosen time steps.

Deciding whether accretion is prograde or retrograde
----------------------------------------------------

We now discuss how to decide whether the sign of spin is positive or negative. In the process of communicating with the inner disk through Lense-Thirring torques, the disk either aligns or counteraligns with the BH spin vector. The condition for which of the two occurs can be derived by assuming that the magnitude of the spin does not change during this alignment (`King et al. 2005 <https://ui.adsabs.harvard.edu/abs/2005MNRAS.363...49K/abstract>`_). Accretion is retrograde if

.. math::
    \cos \theta<-\frac{J_{\mathrm{warp}}}{2 J_{\mathrm{BH}}},
    
where :math:`\cos \theta=\hat{J}_\mathrm{BH}\cdot\hat{J}_\mathrm{d}` is the angle between the initial spin vector and the large-scale angular momentum of the disk. If this condition is not fulfilled, accretion is assumed to be prograde. Note that retrograde accretion is only possible if the angle between the spin vector and the large-scale accretion disk is larger than :math:`90^\circ`, and if the warp angular momentum is comparable to the BH one.

Structure of the warped and precessing accretion disk
-----------------------------------------------------

As mentioned already, Lense-Thirring torques have different effects depending on the type of accretion disk. In particular, their effects depend on the ratio of the viscosity parameter :math:`\alpha` and the aspect ratio :math:`H/R`. For thin discs (:math:`\alpha\gg H/R`), the disk is exactly warped as in the manner described in the preceeding two sections (`Bardeen & Peterson 1975 <https://ui.adsabs.harvard.edu/abs/1975ApJ...195L..65B/abstract>`_). The radius :math:`R_\mathrm{warp}` which separates the inner and outer accretion disc can be calculated by equating the Lense-Thirring precession time-scale (:math:`t_\mathrm{p}=2\pi/\Omega_\mathrm{p}`, with :math:`\Omega_\mathrm{p}=2GJ_\mathrm{BH}/c^2R^3` the precession rate) and the vertical warp propagation time-scale (:math:`t_\mathrm{warp}=R^2/\nu_2`, with :math:`\nu_2` the kinematic viscosity in the vertical direction) (e.g. `Martin et al. 2007 <https://ui.adsabs.harvard.edu/abs/2007MNRAS.381.1617M/abstract>`_). The vertical kinematic viscosity :math:`\nu_2` can be related to the horizontal one, :math:`\nu_1`, by :math:`\nu_2=\xi\nu_1`, with :math:`\xi` a numerical parameter given by

.. math::
    \xi=\frac{2}{\alpha^2}\frac{1+7\alpha^2}{4+\alpha^2}

(`Ogilvie 1999 <https://ui.adsabs.harvard.edu/abs/1999MNRAS.304..557O/abstract>`_, see also `Lodato et al. 2010 <https://ui.adsabs.harvard.edu/abs/2010MNRAS.405.1212L/abstract>`_ for a detailed discussion). We use the relation :math:`\dot{M}=3\pi\nu_1 \Sigma` to calculate :math:`\nu_1`, and therefore :math:`\nu_2`. The warp radius will depend on which region of the thin disc we assume, with each having its own expression for :math:`\Sigma`. In region b) of the `Shakura & Sunyaev (1973) <https://ui.adsabs.harvard.edu/abs/1973A%26A....24..337S/abstract>`_ thin disk, the surface density can be expressed as

.. math::
    \Sigma_\mathrm{TD,b}=6.84 \times 10^{5} \mathrm{~g} \mathrm{~cm}^{-2} \alpha^{-4 / 5} \dot{m}^{3 / 5}\left(\frac{M_{\mathrm{BH}}}{10^{8} M_{\odot}}\right)^{1 / 8}\left(\frac{R}{R_{\mathrm{S}}}\right)^{-3 / 5},
    
while in region c) we have

.. math::
    \Sigma_\mathrm{TD,c}=3.41 \times 10^{4} \mathrm{~g} \mathrm{~cm}^{-2} \alpha^{-4 / 5} \dot{m}^{7/10}\left(\frac{M_{\mathrm{BH}}}{10^{8} M_{\odot}}\right)^{1 / 20}\left(\frac{R}{R_{\mathrm{S}}}\right)^{-3 / 4}.
    
These relations lead to the following expressions for :math:`R_\mathrm{warp}`:

.. math::
    R_{\text {warp,TD,b}}=3410 R_{S} a^{5 / 8} \xi^{-5/8}\alpha^{-1 / 2} \dot{m}^{-1 / 4}\left(\frac{M_{\mathrm{BH}}}{10^{8} M_{\odot}}\right)^{1 / 8}
    
(in region b) and

.. math::
    R_\mathrm{warp,TD,c}=2629R_\mathrm{S}a^{4/7}\xi^{-4/7}\alpha^{-16/35}\dot{m}^{-6/35}\bigg(\frac{M_\mathrm{BH}}{10^8\hspace{0.5mm}\mathrm{M}_\odot}  \bigg)^{4/35},
    
(in region c), with :math:`R_\mathrm{S}=2R_\mathrm{G}` the Schwarzschild radius. These warp radii are generally of order :math:`\approx1000R_\mathrm{G}`, which can lead to fairly quick alignment of the thin disk with the large-scale angular momentum direction (quicker than any significant evolution in mass or spin magnitude, illustrating why the inclusion of the effects of Lense-Thirring torques is important).

In the context of thin disks, there is a futher complication. The self-gravity of the disk may become important at large radii (see `Lodato 2007 <https://www.sif.it/riviste/sif/ncr/econtents/2007/030/07/article/0>`_ for a review). The disk will fragment in the region where the Toomre parameter is :math:`Q(R)>1`. We thus assume that the disk extends out to where :math:`Q(R_\mathrm{sg})=1`. The self-gravity radius :math:`R_\mathrm{sg}` can be calculated from this condition and the definition of the Toomre parameter :math:`Q=\Omega c_{\mathrm{s}} /(\pi G \Sigma)`, yielding

.. math::
    R_{\text {sg,TD,b}}=6460 R_{S} \alpha^{28/51} \dot{m}^{-18/51}\left(\frac{M_{\mathrm{BH}}}{10^{8} M_{\odot}}\right)^{-49/51}
    
in region b) and

.. math::
    R_\mathrm{sg,TD,c}=2456 R_{S} \alpha^{28/45} \dot{m}^{-22/45}\left(\frac{M_{\mathrm{BH}}}{10^{8} M_{\odot}}\right)^{-52/45}
    
in region c). In all our calculations involving :math:`R_\mathrm{warp}` (for deciding the sign of spin and evolving the direction of angular momentum, as described in the preceeding sections), we always take the minimum of :math:`R_\mathrm{warp}` and :math:`R_\mathrm{sg}`. This is because if :math:`R_\mathrm{sg}<R_\mathrm{warp}`, the entire disk of extent :math:`R_\mathrm{sg}` will be warped.

The thick disk does not experience the Bardeen-Peterson effect, i.e. it is never truly aligned nor counter-aligned in its inner regions. Instead, the disk precesses out to several dozen :math:`R_\mathrm{G}`, as seen in simulations (e.g. `Fragile et al. 2007 <https://ui.adsabs.harvard.edu/abs/2007ApJ...668..417F/abstract>`_), and likely observations through quasi-periodic oscillations (QPO; e.g. `Ingram et al. 2012 <https://ui.adsabs.harvard.edu/abs/2012MNRAS.419.2369I/abstract>`_). The slim disk has received much less attention in both simulations and observations (it is both harder to simulate and observe), but its similarity to the thick disk in its geometric aspects likely means that it precesses in a similar manner.

The exact behaviour of the thick and slim disk (which we will collectively call the advection-dominated disks) again depends on the ratio of :math:`\alpha` and :math:`H/R`. Unfortunately, the advection-dominated disks both satisfy :math:`\alpha\approx H/R`, and in this regime, the effects of Lense-Thirring torques are not well understood from a theoretical perspective. However, if :math:`\alpha\ll H/R` (the so-called bending-wave regime), Lense-Thirring torques are known to cause precession of the entire inner disk as a solid body, as seen in observations and simulations. For simplicity, we will thus assume this to be the case for advection-dominated disks.

`Lubow et al. (2002) <https://ui.adsabs.harvard.edu/abs/2002MNRAS.337..706L/abstract>`_ studied the bending-wave regime. In the inner regions, the disk precesses around the spin axis, while in the outer regions, it is aligned with the large-scale angular momentum of the disk. Based on their results the transition radius between the precessing and non-precessing regions of the disk given by

.. math::
    R_\mathrm{warp,adv}=R_\mathrm{G}\bigg(\frac{384a}{25(H/R)^2}\bigg)^{2/5}.
    
In our model, we assume that the inner regions of the disks are on average aligned or counteraligned with the spin vector (one can think of this as averaging over the precession, which has periods of :math:`\approx`days, over long enough time scales). For simplicity, we  also refer to the radii within which this is true as the warp radii. For both of the advection-dominated disks, these radii are only of order several :math:`R_\mathrm{G}`. Note that similar values are found if one assumes that the Bardeen-Peterson effect operates in these disks. While there are some uncertainties in the assumptions we have made, we point out that using any of these values is much more physically motivated than using thin disk equations (the warp radii of order thousands of :math:`R_\mathrm{G}`), which is what is often done (e.g. `Griffin et al. 2019a <https://ui.adsabs.harvard.edu/abs/2019MNRAS.487..198G/abstract>`_, `Dubois et al. 2012 <https://ui.adsabs.harvard.edu/abs/2014MNRAS.440.1590D/abstract>`_).

In order to determine the sign of spin and evolve the angular momentum direction, expressions for the warp mass :math:`M_\mathrm{warp}` and warp angular momentum :math:`J_\mathrm{warp}` are also needed. We calculate this using surface integrals as

.. math::
    M_\mathrm{warp}(R_\mathrm{warp})=2\pi\int_0^{R_\mathrm{warp}}\Sigma(R)R\mathrm{d}R,
    
and
    
.. math::
    J_\mathrm{warp}(R_\mathrm{warp})=2\pi\int_0^{R_\mathrm{warp}}L(R)\Sigma(R)R\mathrm{d}R,
    
respectively. Here, :math:`L(R)` is the specific angular momentum. In the case of the thin disk, we assume Keplerian orbits, i.e. :math:`L(R)=\sqrt{M_\mathrm{BH}G R}`. For the advection-dominated disks, we assume that they are smaller by a numerical factor :math:`\Omega_0`, which is given in the self-similar solutions for the thick (`Narayan & Yi 1995b <https://ui.adsabs.harvard.edu/abs/1995ApJ...452..710N/abstract>`_) and slim disk (`Wang & Zhou 1999 <https://ui.adsabs.harvard.edu/abs/1999ApJ...516..420W/abstract>`_), seperately. The surface densities in both of these accretion disks are given by the same formula in the self-similar solutions, which is

.. math::
    \Sigma_\mathrm{adv}=\frac{\dot{M}}{2\pi R\vert v_\mathrm{r} \vert},

where :math:`v_\mathrm{r}=-\alpha v_0 v_\mathrm{K}` is the radial velocity. Here, :math:`v_\mathrm{K}=\sqrt{M_\mathrm{BH}G/R}` is the Keplerian velocity, and :math:`v_0` is another numerical coefficient that differs between the two solutions. In the thick disk, the numerical coefficients are given by :math:`v_0=3/(5+2\varepsilon)` and :math:`\Omega_0=\sqrt{2\varepsilon/(5+2\varepsilon)}`, where :math:`\varepsilon=(5/3-\gamma)/(\gamma-1)`. The adiabatic index depends on how magnetized the disk is. In particular, it depends on the gas-to-total pressure ratio as :math:`\gamma = (8-3\beta)/(6-3\beta)`, and :math:`\beta` itself depends on :math:`\alpha` (see discussion above on radiative efficiency in the thin disk). :math:`v_0` varies weakly with :math:`\alpha`; for :math:`\alpha=0.05`, it is :math:`0.56`, whereas for :math:`\alpha=0.3`, it evaluates to 0.5. :math:`\Omega_0` depends on :math:`\alpha` somewhat more strongly; we obtain :math:`0.27` and :math:`0.41` for the same values of :math:`\alpha`. The latter value agrees well with the ratio of actual to Keplerian (ISCO) orbital velocity at the event horizon, which is :math:`0.45`. For the slim disc, :math:`v_0=\Omega_0=1/\sqrt{\gamma}`, with :math:`\gamma=5`.

Black hole mergers
------------------

In the process of merging, BHs interact in a very complicated manner. Their final spin is not trivial to predict, and it can depend on a very large parameter space (including the mass ratio of the black holes and the relative orientation and magnitude of the spins). Orbital angular momentum plays a role in the merger as well. We use the fitting function found by `Rezzolla et al. (2007) <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.78.044002>`_, whose results have been found to be very accurate in newer and more sophisticated studies that sweep the huge parameter space of possible merger configurations. The only flaw in these formulas is that they do not include the effects of gravitational radiation. However, the effects of this radiation is confined to a :math:`\approx10\%` level, and only if either of the spin vectors is aligned or counteraligned with the direction of the orbital angular momentum (if it is not, the fits are even more accurate).

The final spin, according to `Rezzolla et al. (2007) <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.78.044002>`_ can be calculated as

.. math::
    \mathbf{a}_\mathrm{fin} = \frac{1}{(1+q)^2}(\mathbf{a}_1+\mathbf{a}_2q^2+\mathbf{l}q),

where :math:`q=M_2/M_1` is the mass ratio (such that :math:`M_2<M_1`), :math:`\mathbf{a}_1` and :math:`\mathbf{a}_2` are the spin vectors, and :math:`\mathbf{l}` is a vector whose direction is the same as that of the orbital angular momentum :math:`\mathbf{L}` (in the centre-of-mass frame), while its magnitude is given by

.. math::
    |\mathbf{l}|=\frac{s_{4}}{\left(1+q^{2}\right)^{2}}\left(\left|\mathbf{a}_{1}\right|^{2}+|\mathbf{a}|_{1}^{2} q^{4}+2\left|\mathbf{a}_{1} \| \mathbf{a}_{2}\right| q^{2} \cos \phi\right)+ \\
    \left(\frac{s_{5} \mu+t_{0}+2}{1+q^{2}}\right)\left(\left|\mathbf{a}_{1}\right| \cos \theta+\left|\mathbf{a}_{2}\right| q^{2} \cos \xi\right)+ \\
    2 \sqrt{3}+t_{2} \mu+t_{3} \mu^{2}.

Here, :math:`\mu=q/(1+q)^2` is the symmetric mass ratio, and :math:`s_4 = -0.129`, :math:`s_5 = -0.384`, :math:`t_0 = -2.686`, :math:`t_2 = -3.454`, :math:`t_3 = 2.353`. The three cosines depend on the angles between the different vectors which play a role in the merger: :math:`\cos \phi=\hat{\mathbf{a}_{1}} \cdot \hat{\mathbf{a}_{\mathbf{2}}}`, :math:`\cos \theta=\hat{\mathbf{a}_{1}} \cdot \hat{\mathbf{l}}`, :math:`\cos \xi=\hat{\mathbf{a}_{2}} \cdot \hat{\mathbf{l}}`.
