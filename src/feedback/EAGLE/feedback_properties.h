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

struct feedback_props {

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

void feedback_props_init(struct feedback_props *fp,
                         const struct phys_const *phys_const,
                         const struct unit_system *us,
                         struct swift_params *params,
                         const struct hydro_props *hydro_props,
                         const struct cosmology *cosmo);
