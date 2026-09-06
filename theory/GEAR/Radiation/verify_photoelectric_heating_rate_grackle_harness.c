/* Standalone Grackle call-through harness for photoelectric_heating=2.
 * Companion to verify_photoelectric_heating_rate.py, which compiles and
 * runs this file; not part of the SWIFT build.
 *
 * Units: density_units=m_H, length_units=1 cm, time_units=1 s,
 * a_value=a_units=1 (non-cosmological). This makes
 * dom = density_units*a^3/m_H = 1 exactly, so a code-unit mass density
 * numerically equals a number density in cm^-3 for pure hydrogen gas
 * (HydrogenFractionByMass forced to 1.0 to avoid He bookkeeping), and
 * coolunit = length_units^2*m_H/time_units^3 = m_H exactly, giving a
 * direct, traceable erg/s/cm^3 <-> code-unit conversion for edot -- this
 * sidesteps hand-deriving Grackle's comoving dom/coolunit bookkeeping.
 *
 * Isolates the photoelectric contribution by calling
 * calculate_cooling_time() twice per Z', identical except
 * photoelectric_heating=2 vs 0: edot = internal_energy_density /
 * cooling_time, and the difference between the two calls is Gamma_PE
 * (up to a negligible cross-term from other channels barely shifting
 * with edot -- checked to be ~1e-49 erg/s/cm^3, i.e. noise, at this
 * script's test point).
 *
 * internal_energy=T_K/temperature_units is not physical T_K by itself
 * (Grackle's T includes a (gamma-1)*mean-molecular-weight factor); each
 * run_case() calibrates it with one calculate_temperature() correction
 * pass so the reported T_on/T_off are the actual achieved temperature,
 * not the naive input guess.
 *
 * igammah=2's gammaha_eff is pinned to 0 above T=2e4 K
 * (cool1d_multi_g.F): irrelevant at this script's T~100K test point, but
 * anyone reusing this harness at higher T, or for igammah=3 (whose
 * pe_eps is T-dependent), needs to know that cutoff exists.
 *
 * Prints one "DATA,Zprime,edot_on,edot_off,T_on,T_off" CSV line per
 * swept Z' for the Python driver to parse; erg/s/cm^3 and K.
 */

#include <grackle.h>
#include <stdio.h>
#include <stdlib.h>

#define MH_CGS 1.67262171e-24 /* g, matches Grackle's own internal mh */

static int run_case(int photoelectric_heating, double G0, double n_H_cgs,
                    double T_K, double Zprime, double *edot_cgs_out,
                    double *T_achieved_K_out) {
  code_units my_units;
  my_units.comoving_coordinates = 0;
  my_units.density_units = MH_CGS;
  my_units.length_units = 1.0;
  my_units.time_units = 1.0;
  my_units.a_units = 1.0;
  my_units.a_value = 1.0;
  set_velocity_units(&my_units);

  chemistry_data *my_chemistry = malloc(sizeof(chemistry_data));
  if (!set_default_chemistry_parameters(my_chemistry)) {
    fprintf(stderr, "set_default_chemistry_parameters failed\n");
    return 0;
  }

  grackle_data->use_grackle = 1;
  grackle_data->with_radiative_cooling = 1;
  grackle_data->primordial_chemistry =
      1;                           /* H, He only, no data-file dependence */
  grackle_data->metal_cooling = 0; /* no Cloudy metal-line cooling */
  grackle_data->UVbackground = 0;  /* no photoionization heating */
  grackle_data->dust_chemistry = 0;
  grackle_data->h2_on_dust = 0;
  grackle_data->photoelectric_heating = photoelectric_heating;
  grackle_data->use_isrf_field = 1; /* match SWIFT: per-cell isrf_habing */
  grackle_data->use_dust_density_field = 0;   /* dust2gas = fgr * metallicity */
  grackle_data->HydrogenFractionByMass = 1.0; /* pure H: skip He bookkeeping */
  grackle_data->Gamma = 5.0 / 3.0;

  if (!initialize_chemistry_data(&my_units)) {
    fprintf(stderr, "initialize_chemistry_data failed\n");
    return 0;
  }

  grackle_field_data my_fields;
  const int field_size = 1;
  my_fields.grid_rank = 3;
  my_fields.grid_dimension = malloc(3 * sizeof(int));
  my_fields.grid_start = malloc(3 * sizeof(int));
  my_fields.grid_end = malloc(3 * sizeof(int));
  my_fields.grid_dx = 0.0;
  for (int i = 0; i < 3; i++) {
    my_fields.grid_dimension[i] = 1;
    my_fields.grid_start[i] = 0;
    my_fields.grid_end[i] = 0;
  }
  my_fields.grid_dimension[0] = field_size;
  my_fields.grid_end[0] = field_size - 1;

  my_fields.density = malloc(field_size * sizeof(gr_float));
  my_fields.internal_energy = malloc(field_size * sizeof(gr_float));
  my_fields.x_velocity = malloc(field_size * sizeof(gr_float));
  my_fields.y_velocity = malloc(field_size * sizeof(gr_float));
  my_fields.z_velocity = malloc(field_size * sizeof(gr_float));
  my_fields.HI_density = malloc(field_size * sizeof(gr_float));
  my_fields.HII_density = malloc(field_size * sizeof(gr_float));
  my_fields.HeI_density = malloc(field_size * sizeof(gr_float));
  my_fields.HeII_density = malloc(field_size * sizeof(gr_float));
  my_fields.HeIII_density = malloc(field_size * sizeof(gr_float));
  my_fields.e_density = malloc(field_size * sizeof(gr_float));
  my_fields.HM_density = malloc(field_size * sizeof(gr_float));
  my_fields.H2I_density = malloc(field_size * sizeof(gr_float));
  my_fields.H2II_density = malloc(field_size * sizeof(gr_float));
  my_fields.DI_density = malloc(field_size * sizeof(gr_float));
  my_fields.DII_density = malloc(field_size * sizeof(gr_float));
  my_fields.HDI_density = malloc(field_size * sizeof(gr_float));
  my_fields.metal_density = malloc(field_size * sizeof(gr_float));
  my_fields.dust_density = NULL; /* unused: use_dust_density_field=0 */
  my_fields.volumetric_heating_rate = malloc(field_size * sizeof(gr_float));
  my_fields.specific_heating_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_HI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_HeI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_HeII_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_H2_dissociation_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_heating_rate = malloc(field_size * sizeof(gr_float));
  my_fields.isrf_habing = malloc(field_size * sizeof(gr_float));

  const double tiny_number = 1.0e-20;
  const double temperature_units = get_temperature_units(&my_units);

  /* n_H_cgs numerically equals the HI+HII code-unit mass density (dom=1,
   * see file header); gas is overwhelmingly neutral at T=100K. */
  for (int i = 0; i < field_size; i++) {
    my_fields.HI_density[i] = n_H_cgs * (1.0 - tiny_number);
    my_fields.HII_density[i] = n_H_cgs * tiny_number;
    my_fields.HeI_density[i] = tiny_number;
    my_fields.HeII_density[i] = tiny_number;
    my_fields.HeIII_density[i] = tiny_number;
    my_fields.e_density[i] = tiny_number;
    my_fields.HM_density[i] = tiny_number;
    my_fields.H2I_density[i] = tiny_number;
    my_fields.H2II_density[i] = tiny_number;
    my_fields.DI_density[i] = tiny_number;
    my_fields.DII_density[i] = tiny_number;
    my_fields.HDI_density[i] = tiny_number;
    my_fields.density[i] = my_fields.HI_density[i] + my_fields.HII_density[i];

    /* metallicity(i) = metal/density/SolarMetalFractionByMass = Zprime */
    my_fields.metal_density[i] =
        Zprime * grackle_data->SolarMetalFractionByMass * my_fields.density[i];

    my_fields.x_velocity[i] = 0.0;
    my_fields.y_velocity[i] = 0.0;
    my_fields.z_velocity[i] = 0.0;
    my_fields.internal_energy[i] = T_K / temperature_units;
    my_fields.volumetric_heating_rate[i] = 0.0;
    my_fields.specific_heating_rate[i] = 0.0;
    my_fields.RT_HI_ionization_rate[i] = 0.0;
    my_fields.RT_HeI_ionization_rate[i] = 0.0;
    my_fields.RT_HeII_ionization_rate[i] = 0.0;
    my_fields.RT_H2_dissociation_rate[i] = 0.0;
    my_fields.RT_heating_rate[i] = 0.0;
    my_fields.isrf_habing[i] = G0;
  }

  /* internal_energy=T_K/temperature_units is NOT physical T_K: Grackle's
   * T also carries a (gamma-1)*mean-molecular-weight factor omitted
   * above. Calibrate with one linear correction pass via
   * calculate_temperature (mu is ~constant with internal_energy at
   * fixed composition, so one pass suffices to floating-point
   * precision) rather than mislabeling the achieved temperature. */
  gr_float *temperature = malloc(field_size * sizeof(gr_float));
  if (!calculate_temperature(&my_units, &my_fields, temperature)) {
    fprintf(stderr, "calculate_temperature failed\n");
    return 0;
  }
  my_fields.internal_energy[0] *= T_K / temperature[0];
  if (!calculate_temperature(&my_units, &my_fields, temperature)) {
    fprintf(stderr, "calculate_temperature failed\n");
    return 0;
  }
  *T_achieved_K_out = temperature[0];

  gr_float *cooling_time = malloc(field_size * sizeof(gr_float));
  if (!calculate_cooling_time(&my_units, &my_fields, cooling_time)) {
    fprintf(stderr, "calculate_cooling_time failed\n");
    return 0;
  }

  /* Fortran's own "energy" is p2d/(gamma-1) = density*internal_energy
   * (ideal gas, code units); edot = energy/cooling_time. */
  const double energy_code =
      my_fields.density[0] * my_fields.internal_energy[0];
  const double edot_code = energy_code / cooling_time[0];

  /* coolunit = length_units^2*m_H/time_units^3 = m_H here (see header). */
  const double coolunit =
      my_units.length_units * my_units.length_units * MH_CGS /
      (my_units.time_units * my_units.time_units * my_units.time_units);
  *edot_cgs_out = edot_code * coolunit;

  free(my_fields.grid_dimension);
  free(my_fields.grid_start);
  free(my_fields.grid_end);
  free(my_fields.density);
  free(my_fields.internal_energy);
  free(my_fields.x_velocity);
  free(my_fields.y_velocity);
  free(my_fields.z_velocity);
  free(my_fields.HI_density);
  free(my_fields.HII_density);
  free(my_fields.HeI_density);
  free(my_fields.HeII_density);
  free(my_fields.HeIII_density);
  free(my_fields.e_density);
  free(my_fields.HM_density);
  free(my_fields.H2I_density);
  free(my_fields.H2II_density);
  free(my_fields.DI_density);
  free(my_fields.DII_density);
  free(my_fields.HDI_density);
  free(my_fields.metal_density);
  free(my_fields.volumetric_heating_rate);
  free(my_fields.specific_heating_rate);
  free(my_fields.RT_HI_ionization_rate);
  free(my_fields.RT_HeI_ionization_rate);
  free(my_fields.RT_HeII_ionization_rate);
  free(my_fields.RT_H2_dissociation_rate);
  free(my_fields.RT_heating_rate);
  free(my_fields.isrf_habing);
  free(cooling_time);
  free(temperature);
  free(my_chemistry);

  return 1;
}

int main(void) {
  grackle_verbose = 0;

  const double n_H_cgs = 100.0; /* cm^-3, total hydrogen */
  const double T_K = 100.0;
  const double G0 = 1.0;
  const double z_prime_sweep[] = {1.0, 0.1, 3.0};
  const int n_sweep = 3;

  for (int zi = 0; zi < n_sweep; zi++) {
    const double Zprime = z_prime_sweep[zi];

    double edot_on, edot_off, T_on, T_off;
    if (!run_case(2, G0, n_H_cgs, T_K, Zprime, &edot_on, &T_on)) return 1;
    if (!run_case(0, G0, n_H_cgs, T_K, Zprime, &edot_off, &T_off)) return 1;

    printf("DATA,%.10e,%.10e,%.10e,%.10e,%.10e\n", Zprime, edot_on, edot_off,
           T_on, T_off);
  }

  return 0;
}
