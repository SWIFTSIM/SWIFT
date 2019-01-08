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
                                 float H, const struct cosmology *restrict cosmo,
				 const struct stars_props *restrict stars_properties,
				 struct xpart *restrict xp) {

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

  /* Add contribution of pj to normalisation of kernel (IMPROVE COMMENT?) */
  si->omega_normalisation_inv += wj / hydro_get_physical_density(pj,cosmo);

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
 * @param d_energy Energy change
 * @param p Particle we're acting on
 * @param cosmo Cosmology struct
 */
static inline void thermal_feedback(float d_energy, struct part * restrict p, 
				    struct xpart * restrict xp,
				    const struct cosmology * restrict cosmo) {

  float u = hydro_get_physical_internal_energy(p, xp, cosmo);
  hydro_set_physical_internal_energy(p, cosmo, u + d_energy);
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
                                  struct spart *restrict si,
                                  struct part *restrict pj, float a, float H,
				  const struct cosmology *restrict cosmo,
				  const struct stars_props *restrict stars_properties,
				  struct xpart *restrict xp) {
  // ALEXEI: GET RID OF a AND H IN SIGNATURE SINCE THESE CAN BE DERIVED FROM COSMO?

  // ALEXEI: THESE CONSTANTS NEED MOVING ELSEWHERE.
  const float total_energy_SNe = 1; // temporary placeholder. actual value 10^51 erg. need to convert to internal units.
  const float units_factor1 = 1.f, units_factor2 = 1.f;

  float wj;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hj_inv = 1.0f / hj;
  const float uj = r * hj_inv;
  kernel_eval(uj, &wj);

  /* Compute weighting for distributing various properties (TODO: better comment?) */
  // ALEXEI: come up with better name for omega_frac?
  float omega_frac = wj/hydro_get_physical_density(pj,cosmo)/si->omega_normalisation_inv;

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

  /* Update mass from SNIa */
  pj->chemistry_data.mass_from_SNIa += si->to_distribute.chemistry_data.mass_from_SNIa * omega_frac;

  /* Update metal mass fraction from SNIa */
  float current_metal_mass_fraction_from_SNIa = pj->chemistry_data.metal_mass_fraction_from_SNIa * current_mass;
  float new_metal_mass_fraction_from_SNIa = current_metal_mass_fraction_from_SNIa + si->to_distribute.chemistry_data.metal_mass_fraction_from_SNIa * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_from_SNIa = new_metal_mass_fraction_from_SNIa/new_mass;

  /* Update mass from SNII */
  pj->chemistry_data.mass_from_SNII += si->to_distribute.chemistry_data.mass_from_SNII * omega_frac;

  /* Update metal mass fraction from SNII */
  float current_metal_mass_fraction_from_SNII = pj->chemistry_data.metal_mass_fraction_from_SNII * current_mass;
  float new_metal_mass_fraction_from_SNII = current_metal_mass_fraction_from_SNII + si->to_distribute.chemistry_data.metal_mass_fraction_from_SNII * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_from_SNII = new_metal_mass_fraction_from_SNII/new_mass;

  /* Update mass from AGB */
  pj->chemistry_data.mass_from_AGB += si->to_distribute.chemistry_data.mass_from_AGB * omega_frac;

  /* Update metal mass fraction from AGB */
  float current_metal_mass_fraction_from_AGB = pj->chemistry_data.metal_mass_fraction_from_AGB * current_mass;
  float new_metal_mass_fraction_from_AGB = current_metal_mass_fraction_from_AGB + si->to_distribute.chemistry_data.metal_mass_fraction_from_AGB * 
    			       si->to_distribute.mass * omega_frac;
  pj->chemistry_data.metal_mass_fraction_from_AGB = new_metal_mass_fraction_from_AGB/new_mass;

  /* Update momentum */
  for (int i = 0; i < 3; i++) {
    // Do we need to calculate relative velocities here?
    pj->v[i] += si->to_distribute.mass * omega_frac * si->v[i];
  }

  /* Energy feedback */
  float d_energy = si->to_distribute.mass * (si->to_distribute.ejecta_specific_thermal_energy 
     + 0.5*(si->v[0]*si->v[0] + si->v[1]*si->v[1] + si->v[2]*si->v[2]) * cosmo->a2_inv);
  if (stars_properties->continuous_heating) {
    // We're doing ONLY continuous heating
    d_energy += si->to_distribute.mass * si->to_distribute.num_SNIa * total_energy_SNe;
  }
  float d_specific_energy = d_energy * omega_frac / current_mass;

  float heating_probability;
  if (!stars_properties->continuous_heating) {
    // We're doing stochastic heating
    heating_probability = units_factor1 * si->to_distribute.num_SNIa *
                          stars_properties->SNIa_energy_fraction /
                          (stars_properties->deltaT_desired * si->to_distribute.ngb_mass);
    // ALEXEI: CHECK UNITS HERE. Eagle does this update in cgs, we should probably keep it in internal units.
    d_specific_energy = stars_properties->deltaT_desired * stars_properties->temp_to_u_factor;
    if (heating_probability >= 1) {
      d_specific_energy = units_factor2 * si->to_distribute.num_SNIa / si->to_distribute.ngb_mass;
      heating_probability = 1; 
    }
  }

  /* pick random number to see if we do stochastic heating */
  unsigned int seed = pj->id;
  if (rand_r(&seed) < heating_probability) {
    // ALEXEI: As above, check units
    thermal_feedback(d_specific_energy, pj, xp, cosmo);
  }

  /* Add in continuous contribution (TEMPORARY COMMENT: from eagle_do_enrich in eagle_enrich.c) */
  thermal_feedback(d_specific_energy, pj, xp, cosmo);

  /* Decrease the mass of star particle (TO CHECK: WHAT ABOUT INTERNAL ENERGY?); */
  si->mass -= si->to_distribute.mass;

}

