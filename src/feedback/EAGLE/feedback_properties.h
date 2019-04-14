/**
 * @brief Stores AGB and SNII yield tables
 */
struct yield_table {
  /* Yield table mass bins */
  double *mass;

  /* Yield table metallicity bins */
  double *metallicity;

  /* Array to store yield table resampled by IMF mass bins */
  double *yield_IMF_resampled;

  /* Array to store yield table being read in */
  double *yield;

  /* Array to store table of ejecta resampled by IMF mass bins */
  double *ejecta_IMF_resampled;

  /* Array to store table of ejecta being read in */
  double *ejecta;

  /* Array to store table of total mass released resampled by IMF mass bins */
  double *total_metals_IMF_resampled;

  /* Array to store table of total mass released being read in */
  double *total_metals;
};

/**
 * @brief Stores tables to determine stellar lifetimes. Used for calculation of
 * IMF
 */
struct lifetime_table {
  /* number of elements, mass, and initial metallicity bins */
  int n_mass;
  int n_z;

  /* table of masses */
  double *mass;

  /* table of metallicities */
  double *metallicity;

  /* table of lifetimes depending on mass an metallicity */
  double **dyingtime;
};

struct feedback_properties {

  /* Flag to switch between continuous and stochastic heating */
  int continuous_heating;

  /* Desired temperature increase due to supernovae */
  float SNe_deltaT_desired;

  /* Conversion factor from temperature to internal energy */
  float temp_to_u_factor;

  /* Energy released by one supernova */
  float total_energy_SNe;

  /* Kinetic energy of SN ejecta per unit mass (check name with Richard)*/
  float ejecta_specific_thermal_energy;

  /* Solar mass */
  float const_solar_mass;

  /* Flag for testing energy injection */
  int const_feedback_energy_testing;

  /* Yield tables for AGB and SNII  */
  struct yield_table yield_AGB;
  struct yield_table yield_SNII;

  /* Array of adjustment factors for SNII  */
  double *typeII_factor;

  /* Arrays of yield tables for SNIa */
  double *yield_SNIa_IMF_resampled;
  double yield_SNIa_total_metals_IMF_resampled;
  double *yields_SNIa;

  /* Parameters to SNIa enrichment model  */
  int SNIa_mode;
  float SNIa_efficiency;
  float SNIa_timescale_Gyr;

  /* Arrays for names of elements being tracked for each enrichment channel */
  char **SNIa_element_names;
  char **SNII_element_names;
  char **AGB_element_names;

  /* Element name string length */
  int element_name_length;

  /* Sizes of dimensions of arrays in yield tables for each enrichment channel
   */
  int SNIa_n_elements;
  int SNII_n_mass;
  int SNII_n_elements;
  int SNII_n_z;
  int AGB_n_mass;
  int AGB_n_elements;
  int AGB_n_z;

  /* log10 of max and min allowable masses for SNII and SNIa in msun */
  float log10_SNII_min_mass_msun;
  float log10_SNII_max_mass_msun;
  float log10_SNIa_max_mass_msun;

  /* Array of mass bins for yield calculations */
  double *yield_mass_bins;

  /* Parameters for IMF  */

  /* IMF model name */
  char IMF_Model[10];

  /* Exponent for IMF if using power law */
  float IMF_Exponent;

  /* Array to store calculated IMF */
  float *imf;

  /* Arrays to store IMF mass bins */
  float *imf_mass_bin;
  float *imf_mass_bin_log10;

  /* Number of IMF mass bins, maximum and minimum bins */
  int n_imf_mass_bins;
  float imf_max_mass_msun;
  float imf_min_mass_msun;
  float log10_imf_min_mass_msun;
  float log10_imf_max_mass_msun;

  /* Table of lifetime values */
  struct lifetime_table lifetimes;

  /* Location of yield tables */
  char yield_table_path[50];

  /* number of type II supernovae per solar mass */
  float num_SNII_per_msun;

  /* wind delay time for SNII */
  float SNII_wind_delay;

  /* Value to set birth time of stars read from ICs */
  float spart_first_init_birth_time;
};

/**
 * @brief Initialize the global properties of the feedback scheme.
 *
 * By default, takes the values provided by the hydro.
 *
 * @param sp The #feedback_properties.
 * @param phys_const The physical constants in the internal unit system.
 * @param us The internal unit system.
 * @param params The parsed parameters.
 * @param p The already read-in properties of the hydro scheme.
 */
INLINE static void feedbakc_props_init(struct feedback_props *fp,
                                       const struct phys_const *phys_const,
                                       const struct unit_system *us,
                                       struct swift_params *params,
                                       const struct hydro_props *p,
                                       const struct cosmology *cosmo) {

  /* Read SNIa timscale */
  fp->SNIa_timescale_Gyr =
      parser_get_param_float(params, "EAGLEFeedback:SNIa_timescale_Gyr");

  /* Read the efficiency of producing SNIa */
  fp->SNIa_efficiency =
      parser_get_param_float(params, "EAGLEFeedback:SNIa_efficiency");

  /* Are we doing continuous heating? */
  fp->continuous_heating =
      parser_get_param_int(params, "EAGLEFeedback:continuous_heating_switch");

  /* Set the delay time before SNII occur */
  const float Gyr_in_cgs = 1e9 * 365 * 24 * 3600;
  fp->SNII_wind_delay =
      parser_get_param_float(params, "EAGLEFeedback:SNII_wind_delay_Gyr") *
      Gyr_in_cgs / units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Read the temperature change to use in stochastic heating */
  fp->SNe_deltaT_desired =
      parser_get_param_float(params, "EAGLEFeedback:SNe_heating_temperature_K");
  fp->SNe_deltaT_desired /=
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  /* Set ejecta thermal energy */
  const float ejecta_velocity =
      1.0e6 / units_cgs_conversion_factor(
                  us, UNIT_CONV_SPEED);  // EAGLE parameter is 10 km/s
  fp->ejecta_specific_thermal_energy = 0.5 * ejecta_velocity * ejecta_velocity;

  /* Energy released by supernova */
  fp->total_energy_SNe =
      1.0e51 / units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* Calculate temperature to internal energy conversion factor */
  fp->temp_to_u_factor =
      phys_const->const_boltzmann_k /
      (p->mu_ionised * (hydro_gamma_minus_one)*phys_const->const_proton_mass);

  /* Read birth time to set all stars in ICs to (defaults to -1 to indicate star
   * present in ICs) */
  fp->spart_first_init_birth_time = parser_get_opt_param_float(
      params, "EAGLEFeedback:birth_time_override", -1);

  /* Copy over solar mass */
  fp->const_solar_mass = phys_const->const_solar_mass;
}
