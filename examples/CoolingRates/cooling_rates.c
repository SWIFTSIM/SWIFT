/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

#if defined(COOLING_EAGLE) && defined(CHEMISTRY_EAGLE) && defined(GADGET2_SPH)
#include "cooling/EAGLE/cooling_rates.h"
#include "cooling/EAGLE/cooling_tables.h"

/* Flag used for printing cooling rate contribution from each
 * element. For testing only. Incremented by 1/(number of elements)
 * until reaches 1 after which point append to files instead of
 * writing new file. */
static float print_cooling_rate_contribution_flag = 0;

/**
 * @brief Wrapper function used to calculate cooling rate and dLambda_du.
 * Writes to file contribution from each element to cooling rate for testing
 * purposes (this function is not used when running SWIFT). Table indices
 * and offsets for redshift, hydrogen number density and helium fraction are
 * passed in so as to compute them only once per particle.
 *
 * @param n_h_i Particle hydrogen number density index
 * @param d_n_h Particle hydrogen number density offset
 * @param He_i Particle helium fraction index
 * @param d_He Particle helium fraction offset
 * @param p Particle structure
 * @param cooling #cooling_function_data structure
 * @param cosmo #cosmology structure
 * @param phys_const #phys_const structure
 * @param abundance_ratio Ratio of element abundance to solar
 */
INLINE static double eagle_print_metal_cooling_rate(
    int n_h_i, float d_n_h, int He_i, float d_He, const struct part *restrict p,
    const struct xpart *restrict xp,
    const struct cooling_function_data *restrict cooling,
    const struct cosmology *restrict cosmo, const struct phys_const *phys_const,
    float *abundance_ratio) {

  /* array to store contributions to cooling rates from each of the
   * elements */
  double *element_lambda;
  element_lambda = malloc((eagle_cooling_N_metal + 2) * sizeof(double));

  /* Get the H and He mass fractions */
  const float XH = p->chemistry_data.metal_mass_fraction[chemistry_element_H];

  /* convert Hydrogen mass fraction in Hydrogen number density */
  const double n_h = hydro_get_physical_density(p, cosmo) * XH /
                     phys_const->const_proton_mass *
                     cooling->number_density_to_cgs;

  /* cooling rate, derivative of cooling rate and internal energy */
  double lambda_net = 0.0;
  double u = hydro_get_physical_internal_energy(p, xp, cosmo) *
             cooling->internal_energy_to_cgs;

  /* Open files for writing contributions to cooling rate. Each element
   * gets its own file.  */
  char output_filename[32];
  FILE **output_file = malloc((eagle_cooling_N_metal + 2) * sizeof(FILE *));

  /* Once this flag reaches 1 we stop overwriting and start appending.  */
  print_cooling_rate_contribution_flag += 1.0 / (eagle_cooling_N_metal + 2);

  /* Loop over each element */
  for (int element = 0; element < eagle_cooling_N_metal + 2; element++) {
    sprintf(output_filename, "%s%d%s", "cooling_element_", element, ".dat");
    if (print_cooling_rate_contribution_flag < 1) {
      /* If this is the first time we're running this function, overwrite the
       * output files */
      output_file[element] = fopen(output_filename, "w");
      print_cooling_rate_contribution_flag += 1.0 / (eagle_cooling_N_metal + 2);
    } else {
      /* append to existing files */
      output_file[element] = fopen(output_filename, "a");
    }
    if (output_file == NULL) {
      error("Error opening file!\n");
    }
  }

  /* calculate cooling rates */
  for (int j = 0; j < eagle_cooling_N_metal + 2; j++) element_lambda[j] = 0.0;
  lambda_net = eagle_metal_cooling_rate(
      log10(u), cosmo->z, n_h, abundance_ratio, n_h_i, d_n_h, He_i, d_He,
      cooling, /*dLambdaNet_du=*/NULL, element_lambda);

  /* write cooling rate contributions to their own files. */
  for (int j = 0; j < eagle_cooling_N_metal + 2; j++) {
    fprintf(output_file[j], "%.5e\n", element_lambda[j]);
  }

  for (int i = 0; i < eagle_cooling_N_metal + 2; i++) fclose(output_file[i]);
  free(output_file);
  free(element_lambda);

  return lambda_net;
}

/**
 * @brief Assign particle density and entropy corresponding to the
 * hydrogen number density and internal energy specified.
 *
 * @param p Particle data structure
 * @param cooling Cooling function data structure
 * @param cosmo Cosmology data structure
 * @param internal_const Physical constants data structure
 * @param nh Hydrogen number density (cgs units)
 * @param u Internal energy (cgs units)
 */
void set_quantities(struct part *restrict p, struct xpart *restrict xp,
                    const struct unit_system *restrict us,
                    const struct cooling_function_data *restrict cooling,
                    const struct cosmology *restrict cosmo,
                    const struct phys_const *restrict internal_const, float nh,
                    double u) {

  double hydrogen_number_density =
      nh * pow(units_cgs_conversion_factor(us, UNIT_CONV_LENGTH), 3);
  p->rho = hydrogen_number_density * internal_const->const_proton_mass /
           p->chemistry_data.metal_mass_fraction[chemistry_element_H];

  float pressure = (u * cosmo->a * cosmo->a) *
                   cooling->internal_energy_from_cgs * p->rho *
                   (hydro_gamma_minus_one);
  p->entropy = pressure * (pow(p->rho, -hydro_gamma));
  xp->entropy_full = p->entropy;
}

/**
 * @brief Produces contributions to cooling rates for different
 * hydrogen number densities, from different metals,
 * tests 1d and 4d table interpolations produce
 * same results for cooling rate, dlambda/du and temperature.
 */
int main(int argc, char **argv) {
  // Declare relevant structs
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct chemistry_global_data chem_data;
  struct part p;
  struct xpart xp;
  struct phys_const internal_const;
  struct cooling_function_data cooling;
  struct cosmology cosmo;
  const char *parametersFileName = "./cooling_rates.yml";

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  const int npts = 250;  // number of values for the internal energy at which
                         // cooling rate is evaluated

  // Set some default values
  float redshift = 0.0, log_10_nh = -1;

  // Read options
  int param;
  while ((param = getopt(argc, argv, "z:d:")) != -1) switch (param) {
      case 'z':
        // read redshift
        redshift = atof(optarg);
        break;
      case 'd':
        // read log10 of hydrogen number density
        log_10_nh = atof(optarg);
        break;
      case '?':
        if (optopt == 'z')
          printf("Option -%c requires an argument.\n", optopt);
        else
          printf("Unknown option character `\\x%x'.\n", optopt);
        error("invalid option(s) to cooling_rates");
    }

  // Read the parameter file
  if (params == NULL) error("Error allocating memory for the parameter file.");
  message("Reading runtime parameters from file '%s'", parametersFileName);
  parser_read_file(parametersFileName, params);

  // Init units
  units_init_from_params(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &internal_const);

  // Init chemistry
  chemistry_init(params, &us, &internal_const, &chem_data);
  chemistry_first_init_part(&internal_const, &us, &cosmo, &chem_data, &p, &xp);
  chemistry_print(&chem_data);

  // Init cosmology
  cosmology_init(params, &us, &internal_const, &cosmo);

  // Set redshift and associated quantities
  const float scale_factor = 1.0 / (1.0 + redshift);
  integertime_t ti_current =
      log(scale_factor / cosmo.a_begin) / cosmo.time_base;
  cosmology_update(&cosmo, &internal_const, ti_current);
  message("Redshift is %f", cosmo.z);

  // Init cooling
  cooling_init(params, &us, &internal_const, &cooling);
  cooling_print(&cooling);
  cooling_update(&cosmo, &cooling);

  // Calculate abundance ratios
  float abundance_ratio[(chemistry_element_count + 2)];
  abundance_ratio_to_solar(&p, &cooling, abundance_ratio);

  // extract mass fractions, calculate table indices and offsets
  float XH = p.chemistry_data.metal_mass_fraction[chemistry_element_H];
  float HeFrac =
      p.chemistry_data.metal_mass_fraction[chemistry_element_He] /
      (XH + p.chemistry_data.metal_mass_fraction[chemistry_element_He]);
  int He_i, n_h_i;
  float d_He, d_n_h;
  get_index_1d(cooling.HeFrac, eagle_cooling_N_He_frac, HeFrac, &He_i, &d_He);

  // Calculate contributions from metals to cooling rate
  // open file
  FILE *output_file = fopen("cooling_output.dat", "w");
  if (output_file == NULL) {
    error("Error opening output file!\n");
  }

  // set hydrogen number density
  const float nh = exp(M_LN10 * log_10_nh);

  /* Initial internal energy */
  double u = 1.0e14;

  // set internal energy to dummy value, will get reset when looping over
  // internal energies
  set_quantities(&p, &xp, &us, &cooling, &cosmo, &internal_const, nh, u);
  float inn_h = hydro_get_physical_density(&p, &cosmo) * XH /
                internal_const.const_proton_mass *
                cooling.number_density_to_cgs;
  get_index_1d(cooling.nH, eagle_cooling_N_density, log10(inn_h), &n_h_i,
               &d_n_h);

  // Loop over internal energy
  for (int j = 0; j < npts; j++) {

    // Update the particle with the new values
    set_quantities(&p, &xp, &us, &cooling, &cosmo, &internal_const, nh,
                   pow(10.0, 10.0 + j * 8.0 / npts));

    // New internal energy
    u = hydro_get_physical_internal_energy(&p, &xp, &cosmo) *
        cooling.internal_energy_to_cgs;

    // calculate cooling rates
    const double temperature = eagle_convert_u_to_temp(
        log10(u), cosmo.z, 0, NULL, n_h_i, He_i, d_n_h, d_He, &cooling);

    const double cooling_du_dt = eagle_print_metal_cooling_rate(
        n_h_i, d_n_h, He_i, d_He, &p, &xp, &cooling, &cosmo, &internal_const,
        abundance_ratio);

    // Dump...
    fprintf(output_file, "%.5e %.5e\n", exp(M_LN10 * temperature),
            cooling_du_dt);
  }
  fclose(output_file);
  message("done cooling rates test");

  /* Clean everything */
  cosmology_clean(&cosmo);
  cooling_clean(&cooling);

  free(params);
  return 0;
}

#else

int main(int argc, char **argv) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  message("This test is only defined for the EAGLE cooling model.");
  return 0;
}
#endif
