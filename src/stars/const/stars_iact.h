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
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_density(float r2, const float *dx, float hi, float hj,
                                 struct spart *restrict si,
                                 const struct part *restrict pj, float a,
                                 float H) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  si->density.wcount += wi;
  si->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

#ifdef DEBUG_INTERACTIONS_STARS
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_STARS)
    si->ids_ngbs_density[si->num_ngb_density] = pj->id;
  ++si->num_ngb_density;
#endif
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
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_stars_feedback(float r2, const float *dx, float hi, float hj,
                                  const struct spart *restrict si,
                                  struct part *restrict pj, float a, float H, float omega_normalisation) {
  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  // ALEXEI: is just plain kernel_eval fine here?
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute weighting for distributing various properties (TODO: better comment?) */
  // ALEXEI: come up with better name for omega_frac?
  float omega_frac = wi/hydro_get_physical_density(pj,cosmo)*omega_normalisation;

  /* Update particle mass */
  float current_mass = hydro_get_mass(pj);
  float new_mass = current_mass + si->to_distribute.mass*omega_frac;
  hydro_set_mass(pj,new_mass);

  // ALEXEI: do we want to use the smoothed mass fraction?
  /* Update total metallicity */
  float current_metal_mass_total = pj->chemistry_data.metal_mass_fraction_total * current_mass;
  float new_metal_mass_total = current_metal_mass_total + si->to_distribute.chemistry_data.metal_mass_fraction_total * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_total = new_metal_mass_total/new_mass;
  
  /* Update mass fraction of each tracked element  */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    float current_metal_mass = pj->chemistry_data.metal_mass_fraction[elem] * current_mass;
    float new_metal_mass = current_metal_mass + si->to_distribute.chemistry_data.metal_mass_fraction[elem] * 
    			   si->to_distribute.mass * omega_frac;
    pj->chemistry_data.metal_mass_fraction[elem] = new_metal_mass/new_mass;
  }

  /* Update iron mass fraction from SNIa  */
  float current_iron_from_SNIa_mass = pj->chemistry_data.iron_mass_fraction_from_SNIa * current_mass;
  float new_iron_from_SNIa_mass = current_iron_from_SNIa_mass + si->to_distribute.chemistry_data.iron_mass_fraction_from_SNIa * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.iron_mass_fraction_from_SNIa = new_iron_from_SNIa_mass/new_mass;

  /* Update mass fraction from SNIa (TODO: MAKE SURE IT REALLY IS A FRACTION) */
  float current_mass_fraction_from_SNIa = pj->chemistry_data.mass_from_SNIa * current_mass;
  float new_mass_fraction_from_SNIa = current_mass_fraction_from_SNIa + si->to_distribute.chemistry_data.mass_from_SNIa * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.mass_from_SNIa = new_mass_fraction_from_SNIa/new_mass;

  /* Update metal mass fraction from SNIa */
  float current_metal_mass_fraction_from_SNIa = pj->chemistry_data.metal_mass_fraction_from_SNIa * current_mass;
  float new_metal_mass_fraction_from_SNIa = current_metal_mass_fraction_from_SNIa + si->to_distribute.chemistry_data.metal_mass_fraction_from_SNIa * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_from_SNIa = new_metal_mass_fraction_from_SNIa/new_mass;

  /* Update mass fraction from SNII (TODO: MAKE SURE IT REALLY IS A FRACTION) */
  float current_mass_fraction_from_SNII = pj->chemistry_data.mass_from_SNII * current_mass;
  float new_mass_fraction_from_SNII = current_mass_fraction_from_SNII + si->to_distribute.chemistry_data.mass_from_SNII * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.mass_from_SNII = new_mass_fraction_from_SNII/new_mass;

  /* Update metal mass fraction from SNII */
  float current_metal_mass_fraction_from_SNII = pj->chemistry_data.metal_mass_fraction_from_SNII * current_mass;
  float new_metal_mass_fraction_from_SNII = current_metal_mass_fraction_from_SNII + si->to_distribute.chemistry_data.metal_mass_fraction_from_SNII * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_from_SNII = new_metal_mass_fraction_from_SNII/new_mass;

  /* Update mass fraction from AGB (TODO: MAKE SURE IT REALLY IS A FRACTION) */
  float current_mass_fraction_from_AGB = pj->chemistry_data.mass_from_AGB * current_mass;
  float new_mass_fraction_from_AGB = current_mass_fraction_from_AGB + si->to_distribute.chemistry_data.mass_from_AGB * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.mass_from_AGB = new_mass_fraction_from_AGB/new_mass;

  /* Update metal mass fraction from AGB */
  float current_metal_mass_fraction_from_AGB = pj->chemistry_data.metal_mass_fraction_from_AGB * current_mass;
  float new_metal_mass_fraction_from_AGB = current_metal_mass_fraction_from_AGB + si->to_distribute.chemistry_data.metal_mass_fraction_from_AGB * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_from_AGB = new_metal_mass_fraction_from_AGB/new_mass;
}
