
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
runner_iact_nonsym_feedback_density(float r2, const float *dx, float hi,
                                    float hj, struct spart *restrict si,
                                    const struct part *restrict pj,
                                    const struct cosmology *restrict cosmo,
                                    const struct feedback_props *feedback_props,
                                    const struct xpart *restrict xp,
                                    integertime_t ti_current) {

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pj);

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Add mass of pj to neighbour mass of si  */
  si->feedback_data.ngb_mass += hydro_get_mass(pj);

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */

  const float rho = hydro_get_comoving_density(pj);
  if (rho != 0)
    si->feedback_data.density_weighted_frac_normalisation_inv += wi / rho;

  /* Compute contribution to the density */
  si->rho_gas += mj * wi;
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
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    float r2, const float *dx, float hi, float hj,
    const struct spart *restrict si, struct part *restrict pj,
    const struct cosmology *restrict cosmo,
    const struct feedback_props *restrict feedback_props,
    struct xpart *restrict xp, integertime_t ti_current) {

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Compute weighting for distributing feedback quantities */
  float density_weighted_frac;
  float rho = hydro_get_comoving_density(pj);
  if (rho * si->feedback_data.density_weighted_frac_normalisation_inv != 0) {
    density_weighted_frac =
        wi / (rho * si->feedback_data.density_weighted_frac_normalisation_inv);
  } else {
    density_weighted_frac = 0.f;
  }

  /* Update particle mass */
  const float current_mass = hydro_get_mass(pj);
  const float new_mass = current_mass + si->feedback_data.to_distribute.mass *
                                            density_weighted_frac;

  hydro_set_mass(pj, new_mass);

  /* Update total metallicity */
  const float current_metal_mass_total =
      pj->chemistry_data.metal_mass_fraction_total * current_mass;
  const float new_metal_mass_total =
      current_metal_mass_total +
      si->feedback_data.to_distribute.total_metal_mass * density_weighted_frac;
  pj->chemistry_data.metal_mass_fraction_total =
      new_metal_mass_total / new_mass;

  /* Update mass fraction of each tracked element  */
  for (int elem = 0; elem < chemistry_element_count; elem++) {
    const float current_metal_mass =
        pj->chemistry_data.metal_mass_fraction[elem] * current_mass;
    const float new_metal_mass =
        current_metal_mass + si->feedback_data.to_distribute.metal_mass[elem] *
                                 density_weighted_frac;
    pj->chemistry_data.metal_mass_fraction[elem] = new_metal_mass / new_mass;
  }

  /* Update iron mass fraction from SNIa  */
  const float current_iron_from_SNIa_mass =
      pj->chemistry_data.iron_mass_fraction_from_SNIa * current_mass;
  const float new_iron_from_SNIa_mass =
      current_iron_from_SNIa_mass +
      si->feedback_data.to_distribute.Fe_mass_from_SNIa * density_weighted_frac;
  pj->chemistry_data.iron_mass_fraction_from_SNIa =
      new_iron_from_SNIa_mass / new_mass;

  /* Update mass fraction from SNIa  */
  const float current_mass_from_SNIa =
      pj->chemistry_data.mass_from_SNIa * current_mass;
  const float new_mass_from_SNIa =
      current_mass_from_SNIa +
      si->feedback_data.to_distribute.mass_from_SNIa * density_weighted_frac;
  pj->chemistry_data.mass_from_SNIa = new_mass_from_SNIa / new_mass;

  /* Update metal mass fraction from SNIa */
  const float current_metal_mass_fraction_from_SNIa =
      pj->chemistry_data.metal_mass_fraction_from_SNIa * current_mass;
  const float new_metal_mass_fraction_from_SNIa =
      current_metal_mass_fraction_from_SNIa +
      si->feedback_data.to_distribute.metal_mass_from_SNIa *
          density_weighted_frac;
  pj->chemistry_data.metal_mass_fraction_from_SNIa =
      new_metal_mass_fraction_from_SNIa / new_mass;

  /* Update mass fraction from SNII  */
  const float current_mass_from_SNII =
      pj->chemistry_data.mass_from_SNII * current_mass;
  const float new_mass_from_SNII =
      current_mass_from_SNII +
      si->feedback_data.to_distribute.mass_from_SNII * density_weighted_frac;
  pj->chemistry_data.mass_from_SNII = new_mass_from_SNII / new_mass;

  /* Update metal mass fraction from SNII */
  const float current_metal_mass_fraction_from_SNII =
      pj->chemistry_data.metal_mass_fraction_from_SNII * current_mass;
  const float new_metal_mass_fraction_from_SNII =
      current_metal_mass_fraction_from_SNII +
      si->feedback_data.to_distribute.metal_mass_from_SNII *
          density_weighted_frac;
  pj->chemistry_data.metal_mass_fraction_from_SNII =
      new_metal_mass_fraction_from_SNII / new_mass;

  /* Update mass fraction from AGB  */
  const float current_mass_from_AGB =
      pj->chemistry_data.mass_from_AGB * current_mass;
  const float new_mass_from_AGB =
      current_mass_from_AGB +
      si->feedback_data.to_distribute.mass_from_AGB * density_weighted_frac;
  pj->chemistry_data.mass_from_AGB = new_mass_from_AGB / new_mass;

  /* Update metal mass fraction from AGB */
  const float current_metal_mass_fraction_from_AGB =
      pj->chemistry_data.metal_mass_fraction_from_AGB * current_mass;
  const float new_metal_mass_fraction_from_AGB =
      current_metal_mass_fraction_from_AGB +
      si->feedback_data.to_distribute.metal_mass_from_AGB *
          density_weighted_frac;
  pj->chemistry_data.metal_mass_fraction_from_AGB =
      new_metal_mass_fraction_from_AGB / new_mass;

  /* Update momentum */
  for (int i = 0; i < 3; i++) {
    pj->v[i] += si->feedback_data.to_distribute.mass * density_weighted_frac *
                (si->v[i] - pj->v[i]);
  }

  /* Energy feedback */
  float u_init = hydro_get_physical_internal_energy(pj, xp, cosmo);
  float heating_probability = -1.f, du = 0.f, d_energy = 0.f;
  d_energy += si->feedback_data.to_distribute.d_energy * density_weighted_frac;

  if (feedback_props->continuous_heating) {
    // We're doing ONLY continuous heating
    d_energy += si->feedback_data.to_distribute.num_SNIa *
                feedback_props->total_energy_SNe * density_weighted_frac *
                si->mass_init;
  } else {
    // We're doing stochastic heating
    heating_probability = si->feedback_data.to_distribute.heating_probability;
    du = feedback_props->SNe_deltaT_desired * feedback_props->temp_to_u_factor;

    if (heating_probability >= 1) {
      du = feedback_props->total_energy_SNe *
           si->feedback_data.to_distribute.num_SNIa /
           si->feedback_data.ngb_mass;
      heating_probability = 1;
    }

    double random_num = random_unit_interval(pj->id, ti_current,
                                             random_number_stellar_feedback);
    if (random_num < heating_probability) {
      message(
          "we did some heating! id %llu star id %llu probability %.5e "
          "random_num %.5e du %.5e "
          "du/ini %.5e",
          pj->id, si->id, heating_probability, random_num, du,
          du / hydro_get_physical_internal_energy(pj, xp, cosmo));
      hydro_set_physical_internal_energy(pj, xp, cosmo, u_init + du);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, u_init + du);
    }
  }

  /* Add contribution from thermal and kinetic energy of ejected material
  (and
   * continuous SNIa feedback) */
  u_init = hydro_get_physical_internal_energy(pj, xp, cosmo);
  du = d_energy / hydro_get_mass(pj);
  hydro_set_physical_internal_energy(pj, xp, cosmo, u_init + du);
  hydro_set_drifted_physical_internal_energy(pj, cosmo, u_init + du);
}
