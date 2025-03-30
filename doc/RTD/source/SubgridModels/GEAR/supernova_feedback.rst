.. Supernova feedback in GEAR model
   Darwin Roduit, 30 March 2025

.. sn_feedback_GEAR_model:

.. _sn_feedback_GEAR_model:


GEAR supernova feedback
=======================
..
   \subsection{Motivation} //find a better name

   Note : Add in the references my master thesis. This way I avoid autoplagiat…

   When a star goes into a supernova, we compute the amount of internal energy, mass and metals the explosion transfers to the star's neighbouring gas particles. We will group all these in the “fluxes” term. 
   We have two models for the distribution of these fluxes and the subgrid modelling of the supernovae : GEAR model and GEAR mechanical model. 

   Note : We may sometimes refer GEAR feedback as GEAR thermal feedback to clearly distinguish it from GEAR mechanical feedback.  

   ## GEAR model

   In the GEAR (thermal)  model, the fluxes are distributed by weighing with the \ac{SPH} kernel  

   $w_{sj} = W_i(\left\| \vect{x}_{sj} \right\|, \, h_s) m_j / \rho_s$ \citep{Revaz_2012_GEAR} for $s$ the star and $j$ the gas. 

   In the GEAR model, we do not inject momentum, only internal energy. Then, internal energy conversion to kinetic energy is left to the hydrodynamic solver, which will compute appropriately the gas density, temperature and velocity. \\
   However, if the cooling radius $R_{\text{cool}}$ of the explosion is unresolved, i.e. the cooling radius is smaller than our simulation resolution, the cooling radiates away the internal energy. 

   To understand why this happens, let us remind the main phases of an SN explosion in a homogeneous medium. We provide a simple picture that is more complicated than the one explained here. See \citep{haid_2016, thornton_1998} for further details. 
   %
   \begin{itemize}
     \item The first stage of the \ac{SN} explosion is the free expansion. In this momentum-conserving regime, the ejecta of the stars sweeps the \ac{ISM}. At the end of this phase, $72 \%$ of the initial SN energy has been converted to thermal energy. 
     \item Once the SN ejecta has swept an \ac{ISM} mass of comparable mass, the blast wave enters the energy-conserving Sedov-Taylor phase. It continues with an adiabatic expansion, performing some $P \dd V$ work on the gas. In this phase, the internal energy is converted into kinetic energy as the ejecta continues sweeping the \ac{ISM} gas. This phase continues until radiative losses become significant after some radius $R_{\text{cool}}$. 
     \item At this point, the blast wave enters the momentum-conserving snowplough phase and forms a thin shell. In this regime, efficient cooling radiates away the internal energy, and thus, the blast wave slows down. 
   \end{itemize}

   Now, we better understand why the internal energy is radiated away. It is a consequence of efficient cooling in the snowplough phase. When this happens, the feedback is unresolved and its energy does not affect the ISM, apart from the mass and metal injection. To circumvent this problem, GEAR thermal feedback implements a fixed delayed cooling mechanism. The cooling of the particles affected by feedback is deactivated during some mega year, ususally \qty{5}{\mega \year} in our simulations. This mechanism allows the internal energy to transform into kinetic energy without immediately being radiated away. However, such an approach poses the question of the time required to prevent gas from cooling in the simulations. 

   ## GEAR mechanical model 

   We implemented two mechanical feedback schemes to better model the blast wave expansion by considering the most critical phases of the SN explosion. These two new models are called GEAR mechanical 1 and GEAR mechanical 2. 
   Our implementations are based on the work of \citet{hopkins_sn_2018} and \citet{hopkins_2024}, which use their feedback in Fire 2 \citep{hopkins_fire-2_2018} and Fire 3 \citep{hopkins_fire-3_2023}. These two implementations differ in their treatment of energy distribution energy and, thus, momentum to distribute. They will allow us to eliminate the delayed cooling and ensure that feedback events are resolved.  

   Note : We only have implemented the flux distribution and the flux injection. The stellar evolution is identical to GEAR model. 

   ### {Neighbour finding}

   When a star $s$ is eligible for feedback, it looks for its neighbour $j$, i.e. gas particles within the star's smoothing length $H_s = \gamma_k h_a$. Such interactions are called non-symmetric since we only look for gas particles the \emph{star sees}. This is the neighbour-finding strategy of the \gear{} thermal feedback and the other \swift{} feedback modules. However, in the \gear{} mechanical feedback, we also consider gas particles $j$ \emph{that see the star} $s$, i.e. the star is within their smoothing length $H_j = \gamma_k h_j$. Such interactions are called symmetric. In summary, we consider all neighbour gas particles $j$ such that $r_{sj} \equiv \norm{ \vect{x}_{s} - \vect{x}_j}< H_s$ or $r_{sj} < H_j$. \\
   Symmetric interactions enable us to consider gas particles that are farther away but still see the star and are eligible for feedback. For instance, in our cosmological simulations, star form in dense disk-like regions $\rho > 1000-10000$ \unit{atoms \per \cm^3}, which have small smoothing lengths, $h \lesssim \qty{1}{\parsec}$. In the vertical plane of the disk, gas particles are not seen by the star particles, but they have smoothing lengths encompassing the disk gas particles and stars. With this symmetric interaction, they are guaranteed to receive the feedback. The symmetric interaction also enables better isotropy in the distribution of fluxes, as we present below. \citet{hopkins_sn_2018} validated the symmetric interaction in their tests and reported that not performing so biases the feedback deposition.

   Finally, we implemented a maximal neighbour search radius, $r_{\text{max}} = \qty{2}{\kilo \parsec}$ by default, to prevent costly neighbour finding and unphysical fluxes deposition, e.g. injecting metals at high distance from the star. 

   ### {Isotropic distribution of fluxes}

   Note : we could reproduce Hopkins conclusion. We provide the setup and the analysis scripts in example XXX. 

   We changed the SPH weighting method to solid angle weighted scheme. This scheme is guaranteed to be be isotropic in the star $s$ reference frame located at $\vect{x}_s$. The scalar weights are based on the solid angles $\Delta \Omega_j$ subtended by the gas particles $j$ and defined by $\omega_j = \Delta \Omega_j / 4 \pi$. To construct the solid angles $\Delta \Omega_j$, we construct a set of faces that enclose the source star $s$ with some convex hull. We assign each face a vector-oriented area $\vect{A}_j$. Then, we can construct the scalar weight, and thus the solid angle perceived by $s$, as:
   %
   \begin{equation}
     \omega_j \equiv \frac{1}{2} \left( 1 - \cfrac{1}{\sqrt{1 + \cfrac{\vect{A}_j \cdot \vect{\hat{x}}_{js}}{\pi \norm{\vect{x}_{js}}^2 }}}    \right) \approx \frac{\Delta \Omega_j}{4 \pi} \, ,
     \label{eq:feedback_scalar_weigth}
   \end{equation}
   %
   where $\vect{\hat{x}}_{js}$ is the unit vector oriented as $\vect{x}_{js}$. This obscure formula \citep{hopkins_sn_2018} interpolates between $1/2$ for $\norm{\vect{x}_{js}}^2 \ll \norm{\vect{A}_j \cdot \vect{\hat{x}}_{js}}$ and $\norm{\vect{A}_j \cdot \vect{\hat{x}}_{js}} / (4 \pi \norm{\vect{x}_{js}}^2) $ for $\norm{\vect{x}_{js}}^2 \gg \norm{\vect{A}_j \cdot \vect{\hat{x}}_{js}}$. 


   \begin{figure}[h]
     \includegraphics[scale=0.9]{figures/feedback_isotropy}
     \caption{Illustration of the isotropic distribution of the fluxes in the \gear{} mechanical feedbacks. The coloured points are gas particles. In purple, we highlight the distribution of the fluxes to the gas particle $j$ within its solid angle. In this simplistic example, the faces $\vect{A}_i$ close exactly. \emph{Source:} Roduit Darwin}
     \label{fig:feedback_isotropy}
   \end{figure}


   We need to construct the convex hull and face vectors $\vect{A}_j$. \citet{hopkins_sn_2018} provide the following formula for \ac{SPH}:
   %
   \begin{equation}
     \vect{A}_j =  \left( \frac{1}{\bar{n}_s^2} \diffp{W(\norm{\vect{x}_{js}}, \, H_s)}{{\norm{\vect{x}_{js}}}} + \frac{1}{\bar{n}_j^2} \diffp{W(\norm{\vect{x}_{js}}, \, H_j)}{{\norm{\vect{x}_{js}}}} \right) \cdot \unitvect{x}_{js} \qquad \bar{n}_s \equiv \sum_j W(\vect{x}_{js}, \, H_s) \,.
     \label{eq:2}
   \end{equation}
   %
   Notice that we have $\sum_j \vect{A}_j = 0$ for an exact closing convex hull. 

   Note : We verified in our simulations how close to it we were; we found $\norm{\sum_j \vect{A}_j} \approx 10^{-2}-10^{-3}$, which is good enough to demonstrate proper face construction. (Maybe to be removed)

   \myFigure{fig:feedback_isotropy} illustrates the isotropic distribution scheme. \\
   However, the scalar weights $\omega_j$ are insufficient to ensure isotropy since we are also dealing with vector quantities such as the momentum $\vect{p}$. We need vector weights $\vect{w}_j$ to ensure isotropy. The derivation of those weights is mathematically involved, and we redirect the interested reader to Hopkins an 2018 paper. Here, we only give the formulas. First, we define $\unitvect{x}_{js}^{\pm}$ the unit vector component in the plus or minus $\alpha = x, \, y, \,z$ directions as:
   %
   \begin{align}
     (\unitvect{x}_{js}^{+})^{\alpha} &\equiv \norm{\vect{x}_{js}}^{-1} \max(\vect{x}_{js}^{\alpha}, \; 0)  \qquad \unitvect{x}_{js}^{+} = ((\unitvect{x}_{js}^{+})^{x}, \; (\unitvect{x}_{js}^{+})^{y}, \: (\unitvect{x}_{js}^{+})^{z}) \\
     (\unitvect{x}_{js}^{-})^{\alpha} &\equiv \norm{\vect{x}_{js}}^{-1} \min(\vect{x}_{js}^{\alpha}, \; 0)  \qquad \unitvect{x}_{js}^{-} = ((\unitvect{x}_{js}^{-})^{x}, \; (\unitvect{x}_{js}^{-})^{y}, \: (\unitvect{x}_{js}^{-})^{z}) \\
     \unitvect{x}_{js} \equiv & \frac{\vect{x}_{js}}{\norm{\vect{x}}_{js}} = \sum_{+, \, -} \unitvect{x}_{js}^{\pm} \; .
   \end{align}
   %
   Then, we define $(f_{\pm}^{\alpha})_s$ the star's vector isotropy correction factor in the plus or minus direction:
   %
   \begin{equation}
     (f_{\pm}^{\alpha})_s \equiv \left\{ \frac{1}{2} \left[ 1 + \left( \frac{\displaystyle \sum_k \omega_k \abs{(\unitvect{x}_{ks}^{\mp})^{\alpha}} }{\displaystyle \sum_k \omega_k \abs{(\unitvect{x}_{ks}^{\pm})^{\alpha}}}    \right)^2                \right]             \right\}^{1/2} \; .
   \end{equation}
   %
   The vector weigths $\vect{w}_j$ and the normalized vector weights $\vect{\bar{w}}_j$ are thus defined as:
   %
   \begin{align}
     w_j^{\alpha} &\equiv \omega_j \sum_{+, \, -} (\unitvect{x}_{js}^{\pm})^{\alpha} \, (f_{\pm}^{\alpha})_s \\
     \bar{w}_j^{\alpha} &\equiv \frac{w_j^{\alpha}}{\displaystyle \sum_k \norm{\vect{w}_k}} \; .
   \end{align}
   %
   Those expressions are evaluated in two new neighbour loops. Physically, the normalized vector weigths $\vect{\bar{w}}_j$ account for the asymmetries about the vector $\unitvect{x}_{js}$ in the faces $\vect{A}_j$. Those complex mathematical expressions have the following properties:
   %
   \begin{itemize}
     \item The distribution of the fluxes is isotropic. 
     \item They ensure machine-accurate conservation of the fluxes to be distributed.
     \item The fractional error $\norm{\sum_j \vect{p}_j}/ p_{\text{ej}}$, with $p_{\text{ej}}$ the momentum ejected by the supernova, is independent of the spatial distribution of the neighbours in the kernel. 
   \end{itemize}
   %
   \citet{hopkins_sn_2018} provides a complete discussion of those properties and detailed comparisons within the Fire-1 and Fire-2 simulations. \\



   To read : ——>



   Consider now that the supernova explosion must distribute scalar fluxes $X_{\text{ej}}$ such as the mass $m_{\text{ej}}$, the metals $m_{Z, \text{ej}}$, the total energy $E_{\text{ej}}$, as well as vector norm fluxes $Y_{\text{ej}}$ such as the momentum $p_{\text{ej}}$. The fluxes distributed to the gas particles are defined as:
   %
   \begin{align}
     &\Delta X_j = \norm{\vect{\bar{w}}_j} X_{\text{ej}} \\
     &\Delta \vect{Y}_j = \vect{\bar{w}}_j Y_{\text{ej}} \; .
   \end{align}
   %
   The machine-accurate conservation means:
   %
   \begin{align}
     &\sum \Delta X_j = X_{\text{ej}} \\
     &\sum \norm{\Delta \vect{Y}_j} = Y_{\text{ej}} \\
     &\sum \Delta \vect{Y}_j = \vect{0} \; .
   \end{align}
   %

   Until now, we implicitly worked in the reference frame of the star, i.e. $\vect{x}_s = \vect{0}$, $\vect{v}_s \equiv \diff{\vect{x}_s}{{t}} = \vect{0}$. The distribution of flux is isotropic in the reference frame of the stars. However, we need to consider the star motion and thus boost the fluxes in the laboratory frame to obtain the fluxes $\Delta X_j'$. For the mass, metals and momentum, this is trivial:
   %
   \begin{align}
     \Delta m_j' &\equiv \Delta m_j = \norm{\vect{\bar{w}}_j} m_{\text{ej}} \;, \quad \Delta m_{Z, j}' \equiv \Delta m_{Z, j} = \norm{\vect{\bar{w}}_j} m_{Z, \text{ej}} \\
     \Delta \vect{p}_{js} ' &\equiv \Delta \vect{p}_{js} + \Delta m_j \vect{v}_s
   \end{align}
   %
   For the energy, this depends on the implementation. The main differences are that we ignore the star-gas motion in \gear{} mechanical 1, while in \gear{} mechanical 2, we consider this motion. This also changes $p_{\text{ej}}$.\\
   In both implementation, we verify whether we resolve the Sedov-Taylor phase and inject the correct energy, internal energy and mometum into the surrounding gas. The algorithm depends on whether we include the star-gas motion and thus depends on the implementations. For the following, we write:
   %
   \begin{equation}
       \Delta \vect{p}_{js} \equiv \vect{\bar{w}}_{j} p_{0, s} \; ,
   \end{equation}
   % 
   where $ p_{0, s}$ depends on the particular treatment of the star-gas motion and is not simply $p_{\text{ej}}$. 

   \subsubsection{\gear{} mechanical 1}

   In \gear{} mechanical 1, we have the following fluxes to distribute: $m_{\text{ej}}$, $m_{Z, \text{ej}}$ and $E_{\text{ej}}$. The momentum flux is $p_{\text{ej}} = \sqrt{2 m_{\text{ej}} E_{\text{ej}}}$. The fluxes are given to the gas particle $j$ as:
   %
   \begin{align}
     m_j^{\text{new}} &= m_j +  \Delta m_j' = m_j + \norm{\vect{\bar{w}}_j} m_{\text{ej}} \label{eq:gear_m_1_flux_m}\\
     m_{Z, j}^{\text{new}} &= m_{Z, j} +  \Delta m_{Z, j}' = m_{Z, j} + \norm{\vect{\bar{w}}_j} m_{\text{ej}} \label{eq:gear_m_1_flux_Z} \\
     E_j^{\text{new}} &= E_{\text{kin}}^{\text{new}} +  U_{\text{int}}^{\text{new}} =  E_{\text{kin}} +  U_{\text{int}} + \norm{\vect{\bar{w}}_j} E_{\text{ej}} + \frac{1}{2 \Delta m_j} ( \norm{\Delta \vect{p}_{js}'}^2 - \norm{\Delta \vect{p}_{js}}^2 ) \label{eq:gear_m_1_flux_E}\\ 
     U_{\text{int}}^{\text{new}} &= U_{\text{int}} + \Delta U \, , \quad \Delta U = (E_j^{\text{new}} - E_{\text{kin}}^{\text{new}}) - U_{\text{int}} \label{eq:gear_m_1_flux_U} \\
     \vect{p}_j^{\text{new}} &= \vect{p}_j + \Delta m_j \vect{v}_s +  \vect{\bar{w}}_{j} p_{0, s} \label{eq:gear_m_1_flux_p} \; .
   \end{align}
   %
   Now, we need to define $p_{0, s}$. In high-density regions and/or in low-resolution simulations, we may not be able to resolve the Sedov-Taylor expansion phase. As we explained above, during the latter, the blastwave sweeps the gas and thus does some $P \dd V$ work on the gas. This work converts energy into momentum until reaching the end of the phase, when the cooling becomes efficient at some cooling radius $R_{\text{cool}}$. If we do not resolve the Taylor-Sedov phase, we may give an incorrect amount of momentum and energy into the \ac{ISM}. At the beginning of the snowplough phase, the momentum of a supernova reaches some terminal value $p_t$. It can be written as:
   %
   \begin{equation}
       p_t = p_{t, 0} \; \mathcal{F}_{E}(E) \mathcal{F}_{n}(n)  \mathcal{F}_{Z}(Z)  \mathcal{F}_{\vect{v}} (\vect{v}) \;, \label{eq:p_terminal}
   \end{equation}
   %
   where $\mathcal{F}_{k}$ are functions depending on the total \ac{SN}-frame ejecta energy $E$, the gas density $n$, metallicity $Z$ and velocity field $\vect{v}$. We use the same parametrisation than Fire-3 \citep{hopkins_fire-3_2023}, i.e.
   %
   \begin{align}
     p_{t, 0} &= \qty{200}{\Msun \km \per \second}, \qquad \mathcal{F}_E = \frac{E}{10^{51} \unit{erg}} \; , \qquad  \mathcal{F}_{\vect{v}} = 1 \\
     \mathcal{F}_n &= 2.63, \; \tilde{n} < 10^{-3} \text{ and } \mathcal{F}_n = \tilde{n}^{-0.143}, \; \tilde{n} \geq 10^{-3}, \qquad \tilde{n} \equiv  \frac{n}{\unit{\per \cm^3}} \label{eq:terminal_momentum_density_dependence} \\
     \mathcal{F}_Z &= 2, \; \text{ if } \tilde{z} < 10^{-2} \text{ , } \mathcal{F}_Z = \tilde{z}^{-0.18}, \; \text{ if } 10^{-2}\leq  \tilde{z} \leq 1 \text{ and } \mathcal{F}_Z = \tilde{z}^{-0.12}, \; \text{ if } \tilde{z} > 1,  \quad \tilde{z} \equiv  \frac{Z}{Z_{\odot}}
   \end{align}
   %
   Also, we use $\mathcal{F}_{\vect{v}} = 1$ advised by \citet{hopkins_2024}. \\
   To account for the potentially unresolved Taylor-Sedov phase, we first calculate the momentum that would be coupled to the gas particle if the blastwave were energy-conserving throughout this single element. 



   This momentum can be computed by equating the \ac{SN} energy injected right after the explosion into the gas particle $j$ and the energy of $j$ at the end of the Taylor-Sedov phase:
   %
   \begin{align*}
     & \Delta E_j^{\text{initial}} = \Delta E_j^{\text{End}}  \\
     & \Leftrightarrow \frac{\Delta p_j^2}{2 \Delta m_j} = \frac{p_{j, \text{final}}^2}{2 (m_j + \Delta m_j)} \Leftrightarrow (m_j + \Delta m_j) \Delta p_j^2 = p_{j, \text{final}}^2 \Delta m_j \\
     &\Leftrightarrow  p_{j, \text{final}}^2 = \left( 1 + \frac{m_j}{\Delta m_j} \right) \Delta p_j^2 \\
     & \Rightarrow p_{j, \text{final}} = \sqrt{1 + \frac{m_j}{\Delta m_j}} \Delta p_j \, , \; \Delta p_j =  \norm{\vect{\bar{w}}_j} p_{\text{ej}} \; .
   \end{align*}
   %
   Notice that at the end of the energy conserving phase, the final mass of the gas particle is its initial mass $m_j$ plus the ejected mass $\Delta m_j$. \\
   Then, we compare this momentum to the terminal momentum $p_t$ and assign the momentum to be:
   %
   \begin{equation}
     p_{0,s} = p_{\text{ej}} \min \left(\sqrt{1 + \frac{m_j}{\Delta m_j}}, \; \frac{p_t}{p_{\text{ej}}} \right) \; .
   \end{equation}
   %
   The last thing to do is to couple the correct internal energy when the cooling radius $R_{\text{cool}}$ is unresolved. The cooling radius is determined by the value of $p_t$ since, at the end of the Sedov-Taylor phase, we have $R_{\text{cool}} = R_{\text{Shock, SN}}$ and by conservation of energy, $E_{\text{ej}} = p_{\text{ej}}^2 / (2 m_{\text{ej}}) = p_t^2 / (2 (m_{\text{ej}} + m_{\text{swept}}(R_{\text{cool}}))) = E_{\text{End}}$. Then, using $m_{\text{swept}}(R_{\text{cool}}) = \frac{4}{3} \pi R_{\text{cool}}^3 \rho_{\text{swept}}$ and rearranging the terms, we find:
   %
   \begin{equation}
       R_{\text{cool}} = \left( \frac{3 m_{\text{ej}}}{4 \pi \rho} \right)^{1/3} \left(\frac{p_t^2}{p_{\text{ej}}^2} - 1 \right)^{1/3} \; ,
   \label{eq:cooling_radius}
   \end{equation}
   %
   where $\rho$ is the density. \\
   As the internal energy outside $R_{\text{cool}}$ decays $\propto (r/R_{\text{cool}})^{-6.5}$ \citep{thornton_1998}, if $r_j \equiv \norm{\vect{x}_{js}} > R_{\text{cool}}$, we reduce the internal energy as $\Delta U \leftarrow \Delta U (r_j/R_{\text{cool}})^{-6.5}$. Otherwise, we leave $\Delta U$ unchanged. This concludes this part on the implementation of \gear{} mechanical 1 feedback. \\
   Notice the following: we found a simple analytic formula \ref{eq:cooling_radius} for the cooling radius \footnote{Notice that throughout this work, we use the symbols $\rho$ and $n$ for the density interchangeably. }. This formula provide the following density dependance: $R_{\text{cool}} \propto \rho^{-1/3} p_t^{2/3}$. If we use our relation \eqref{eq:terminal_momentum_density_dependence} for the terminal momentum, we have:
   %
   \begin{equation}
	R_{\text{cool}} \propto \rho^{-1/3} \text{ if } \rho < \qty{e-3}{\per \cm^3} \quad \text{and} \quad  R_{\text{cool}} \propto \rho^{-0.428} \text{ if } \rho \geq \qty{e-3}{\per \cm^3} \;.
   \end{equation}
   %

