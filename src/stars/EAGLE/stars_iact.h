/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param xp Extra particle data
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_density(float r2, const float *dx, float hi, float hj,
                                 struct spart *restrict si,
                                 const struct part *restrict pj, 
				 const struct cosmology *restrict cosmo,
				 const struct stars_props *restrict stars_properties,
				 struct xpart *restrict xp, integertime_t ti_current) {

  /* Get the gas mass. */
  const float mj = pj->mass;

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  
  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Compute contribution to the number of neighbours */
  si->density.wcount += wi;
  si->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Add mass of pj to neighbour mass of si  */
  si->ngb_mass += hydro_get_mass(pj);

  /* Add contribution of pj to normalisation of kernel (TODO: IMPROVE COMMENT?) */
  si->omega_normalisation_inv += wj / hydro_get_physical_density(pj,cosmo);
  
  /* Compute contribution to the density */
  si->rho_gas += mj * wi;

#ifdef DEBUG_INTERACTIONS_STARS
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_STARS)
    si->ids_ngbs_density[si->num_ngb_density] = pj->id;
  ++si->num_ngb_density;
#endif
}

/**
 * @brief Increases thermal energy of particle due 
 * to feedback by specified amount
 *
 * @param du change in internal energy
 * @param p Particle we're acting on
 * @param xp Extra particle data
 * @param cosmo Cosmology struct
 */
static inline void thermal_feedback(float du, struct part * restrict p, 
				    struct xpart * restrict xp,
				    const struct cosmology * restrict cosmo) {

  float u = hydro_get_physical_internal_energy(p, xp, cosmo);
  hydro_set_physical_internal_energy(p, cosmo, u + du);
  // Just setting p->entropy is not enough because xp->entropy_full gets updated with p->entropy_dt
  // TODO: ADD HYDRO FUNCTIONS FOR UPDATING DRIFTED AND NON DRIFTED INTERNAL ENERGY AND GET RID OF 
  // THE ENTROPY UPDATE HERE.
  xp->entropy_full = p->entropy;
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle.
 * @param pj Second (gas) particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param xp Extra particle data
 * @param ti_current Current integer time used value for seeding random number generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_feedback(float r2, const float *dx, float hi, float hj,
                                  struct spart *restrict si,
                                  struct part *restrict pj,
				  const struct cosmology *restrict cosmo,
				  const struct stars_props *restrict stars_properties,
				  struct xpart *restrict xp, integertime_t ti_current) {
  float wj;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hj_inv = 1.0f / hj;
  const float uj = r * hj_inv;
  kernel_eval(uj, &wj);

  /* Compute weighting for distributing feedback quantities */
  const float omega_frac = wj/(hydro_get_physical_density(pj,cosmo) * si->omega_normalisation_inv);

  /* Update particle mass */
  const float current_mass = hydro_get_mass(pj);
  float new_mass = current_mass + si->to_distribute.mass*omega_frac;
  if (stars_properties->const_feedback_energy_testing) new_mass = current_mass;
  hydro_set_mass(pj,new_mass);
  
  /* Decrease the mass of star particle */
  si->mass -= si->to_distribute.mass*omega_frac;

  // ALEXEI: do we want to use the smoothed mass fraction?
  /* Update total metallicity */
  const float current_metal_mass_total = pj->chemistry_data.metal_mass_fraction_total * current_mass;
  const float new_metal_mass_total = current_metal_mass_total + si->to_distribute.chemistry_data.metal_mass_fraction_total * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_total = new_metal_mass_total/new_mass;
  
  /* Update mass fraction of each tracked element  */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const float current_metal_mass = pj->chemistry_data.metal_mass_fraction[elem] * current_mass;
    const float new_metal_mass = current_metal_mass + si->to_distribute.chemistry_data.metal_mass_fraction[elem] * 
    			   si->to_distribute.mass * omega_frac;
    pj->chemistry_data.metal_mass_fraction[elem] = new_metal_mass/new_mass;
  }

  /* Update iron mass fraction from SNIa  */
  const float current_iron_from_SNIa_mass = pj->chemistry_data.iron_mass_fraction_from_SNIa * current_mass;
  const float new_iron_from_SNIa_mass = current_iron_from_SNIa_mass + si->to_distribute.chemistry_data.iron_mass_fraction_from_SNIa * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.iron_mass_fraction_from_SNIa = new_iron_from_SNIa_mass/new_mass;

  /* Update mass from SNIa */
  pj->chemistry_data.mass_from_SNIa += si->to_distribute.chemistry_data.mass_from_SNIa * omega_frac;

  /* Update metal mass fraction from SNIa */
  const float current_metal_mass_fraction_from_SNIa = pj->chemistry_data.metal_mass_fraction_from_SNIa * current_mass;
  const float new_metal_mass_fraction_from_SNIa = current_metal_mass_fraction_from_SNIa + si->to_distribute.chemistry_data.metal_mass_fraction_from_SNIa * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_from_SNIa = new_metal_mass_fraction_from_SNIa/new_mass;

  /* Update mass from SNII */
  pj->chemistry_data.mass_from_SNII += si->to_distribute.chemistry_data.mass_from_SNII * omega_frac;

  /* Update metal mass fraction from SNII */
  const float current_metal_mass_fraction_from_SNII = pj->chemistry_data.metal_mass_fraction_from_SNII * current_mass;
  const float new_metal_mass_fraction_from_SNII = current_metal_mass_fraction_from_SNII + si->to_distribute.chemistry_data.metal_mass_fraction_from_SNII * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_from_SNII = new_metal_mass_fraction_from_SNII/new_mass;

  /* Update mass from AGB */
  pj->chemistry_data.mass_from_AGB += si->to_distribute.chemistry_data.mass_from_AGB * omega_frac;

  /* Update metal mass fraction from AGB */
  const float current_metal_mass_fraction_from_AGB = pj->chemistry_data.metal_mass_fraction_from_AGB * current_mass;
  const float new_metal_mass_fraction_from_AGB = current_metal_mass_fraction_from_AGB + si->to_distribute.chemistry_data.metal_mass_fraction_from_AGB * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_from_AGB = new_metal_mass_fraction_from_AGB/new_mass;

  /* Update momentum */
  for (int i = 0; i < 3; i++) {
    // Do we need to calculate relative velocities here?
    pj->v[i] += si->to_distribute.mass * omega_frac * si->v[i];
  }

  /* Energy feedback */
  float heating_probability = -1.f, du = 0.f, d_energy = 0.f;
  d_energy = si->to_distribute.mass * (si->to_distribute.ejecta_specific_thermal_energy 
     + 0.5*(si->v[0]*si->v[0] + si->v[1]*si->v[1] + si->v[2]*si->v[2]) * cosmo->a2_inv);

  // If statement temporary for testing, in practice would always be on.
  if (stars_properties->const_feedback_energy_testing) {
    if (stars_properties->continuous_heating) {
      // We're doing ONLY continuous heating 
      d_energy += si->to_distribute.num_SNIa * stars_properties->total_energy_SNe * omega_frac * si->mass_init;
      du = d_energy/hydro_get_mass(pj);
      thermal_feedback(du,pj,xp,cosmo);
    } else {
      // We're doing stochastic heating
      heating_probability = stars_properties->SNe_temperature_h * si->to_distribute.num_SNIa *
                            stars_properties->SNIa_energy_fraction /
                            (stars_properties->SNe_deltaT_desired * si->ngb_mass);
      du = stars_properties->SNe_deltaT_desired * stars_properties->temp_to_u_factor;
      if (heating_probability >= 1) {
        du = stars_properties->SNe_energy_h * si->to_distribute.num_SNIa / si->ngb_mass;
        heating_probability = 1; 
      }
    }
  }

  /* pick random number to see if we do stochastic heating */
  // Temporary assignment of random seed. Discuss with Matthieu for better 
  // way of generating random numbers
  unsigned int seed = 3*(pj->id + ti_current) % 8191;
  double random_num = rand_r(&seed) * stars_properties->inv_rand_max;
  if (random_num < heating_probability) {
    message("we did some heating! id %llu probability %.5e random_num %.5e du %.5e du/ini %.5e", pj->id, heating_probability, random_num, du, du/hydro_get_physical_internal_energy(pj,xp,cosmo));
    thermal_feedback(du, pj, xp, cosmo);
  }

}
