/*******************************************************************************
 * This file is part of SWIFT.
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

/* Key arrival harness of the geajay cooling module: read a parameter file with
 * the module's own cooling_init() and print what reaches the geajay library.
 * It is run by testGeajayKeys.sh, one parameter file at a time.
 *
 * The library entry points that take options are interposed: the executable
 * defines them, records the arguments and calls the real function.
 * Usage: testGeajayKeys <parameter file> or testGeajayKeys --defaults. */
#define _GNU_SOURCE

/* Config parameters. */
#include <config.h>

#include <dlfcn.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "swift.h"

#include "geajay_cosmic_ray.h"
#include "geajay_options.h"

/* ---- dump of an options struct ---------------------------------------- */
#define P_INT(n) printf("LIB %s=%d\n", #n, o->n);
#define P_DOUBLE(n) printf("LIB %s=%.10g\n", #n, o->n);
#define P_STRING(n) printf("LIB %s=%s\n", #n, o->n ? o->n : "(null)");
#define GEAJAY_OPTION(n, g, t) P_##t(n)
#define GEAJAY_BACKGROUND_OPTION(n, t, d) P_##t(n)

static void dump_options(const char *via, const geajay_options *o) {
  printf("LIBPATH %s\n", via);
#include "geajay_options_fields.def"
#include "geajay_background_fields.def"
}

static void dump_chemistry(const geajay_chemistry_data *c) {
  geajay_options o;
  geajay_options_set_defaults(&o);
  if (geajay_options_from_chemistry_data(c, &o) == GEAJAY_FAIL)
    printf("LIBPATH chemistry_data (options_from_chemistry_data FAILED)\n");
  else
    dump_options("geajay_context_create(chemistry_data)", &o);
}

/* ---- interposed entry points ------------------------------------------ */
geajay_context *geajay_context_create(
    const geajay_chemistry_data *chemistry,
    const geajay_chemistry_data_storage *rates, const geajay_code_units *units,
    geajay_context_error *error) {
  static geajay_context *(*real)(const geajay_chemistry_data *,
                                 const geajay_chemistry_data_storage *,
                                 const geajay_code_units *,
                                 geajay_context_error *);
  if (!real) real = dlsym(RTLD_NEXT, "geajay_context_create");
  dump_chemistry(chemistry);
  return real(chemistry, rates, units, error);
}

geajay_context *geajay_context_create_from_options(
    const struct geajay_options *options, const geajay_code_units *units,
    geajay_context_error *error) {
  static geajay_context *(*real)(const struct geajay_options *,
                                 const geajay_code_units *,
                                 geajay_context_error *);
  if (!real) real = dlsym(RTLD_NEXT, "geajay_context_create_from_options");
  dump_options("geajay_context_create_from_options", options);
  return real(options, units, error);
}

int geajay_context_set_redshift(geajay_context *context, double redshift) {
  static int (*real)(geajay_context *, double);
  if (!real) real = dlsym(RTLD_NEXT, "geajay_context_set_redshift");
  printf("LIBCALL set_redshift=%.10g\n", redshift);
  return real(context, redshift);
}

void geajay_context_set_subcycle_on_failure(geajay_context *context, int n) {
  static void (*real)(geajay_context *, int);
  if (!real) real = dlsym(RTLD_NEXT, "geajay_context_set_subcycle_on_failure");
  printf("LIBCALL set_subcycle_on_failure=%d\n", n);
  real(context, n);
}

int geajay_set_table_density_edge(geajay_chemistry_data_storage *rates,
                                  geajay_table_edge edge) {
  static int (*real)(geajay_chemistry_data_storage *, geajay_table_edge);
  if (!real) real = dlsym(RTLD_NEXT, "geajay_set_table_density_edge");
  printf("LIBCALL set_table_density_edge=%d\n", (int)edge);
  return real(rates, edge);
}

geajay_cosmic_ray_reduced_context *geajay_cosmic_ray_reduced_context_create(
    const geajay_context *context, const geajay_cosmic_ray_settings *s,
    geajay_context_error *error) {
  static geajay_cosmic_ray_reduced_context *(*real)(
      const geajay_context *, const geajay_cosmic_ray_settings *,
      geajay_context_error *);
  if (!real)
    real = dlsym(RTLD_NEXT, "geajay_cosmic_ray_reduced_context_create");
  printf("LIBCALL cosmic_ray_reduced_context_create\n");
  printf("LIB cr.network=%d\nLIB cr.input=%d\nLIB cr.heating=%d\n", s->network,
         s->input, s->heating);
  printf("LIB cr.rate_per_energy_density=%.10g\n", s->rate_per_energy_density);
  printf("LIB cr.constant_rate=%.10g\n", s->constant_rate);
  printf("LIB cr.heat_per_ionisation=%.10g\n", s->heat_per_ionisation);
  printf("LIB cr.coulomb_loss_coefficient=%.10g\n", s->coulomb_loss_coefficient);
  printf("LIB cr.hadronic_loss_coefficient=%.10g\n",
         s->hadronic_loss_coefficient);
  return real(context, s, error);
}

/* ---- module state ------------------------------------------------------ */
#define M_INT(n) printf("MOD %s=%d\n", #n, (int)c.n);
#define M_DBL(n) printf("MOD %s=%.10g\n", #n, (double)c.n);
#define M_STR(n) printf("MOD %s=%s\n", #n, c.n);

int main(int argc, char **argv) {
  if (argc < 2) return 2;
  if (strcmp(argv[1], "--defaults") == 0) {
    geajay_options o;
    geajay_options_set_defaults(&o);
    dump_options("library defaults", &o);
    return 0;
  }
  struct swift_params *params = malloc(sizeof(struct swift_params));
  struct unit_system us;
  struct phys_const phys_const;
  struct hydro_props hydro_properties;
  struct cooling_function_data c;
  memset(&c, 0, sizeof(c));

  parser_read_file(argv[1], params);
  units_init_from_params(&us, params, "InternalUnitSystem");
  phys_const_init(&us, params, &phys_const);
  hydro_props_init(&hydro_properties, &phys_const, &us, params);

  cooling_init(params, &us, &phys_const, &hydro_properties, &c);
  fflush(stdout);

  /* The chemistry_data after geajay_initialize_chemistry_data(), as module
     fields. The options of the context are above. */
  printf("LIBPATH done\n");
  M_INT(primordial_chemistry) M_INT(with_uv_background)
  M_INT(with_metal_cooling) M_INT(H2_three_body_rate) M_INT(H2_cie_cooling)
  M_INT(dust_model) M_INT(dust_to_gas_model) M_INT(photoelectric_heating_model)
  M_DBL(photoelectric_heating_rate) M_INT(dust_recombination_cooling)
  M_INT(H2_on_dust) M_INT(H2_photodissociation_heating)
  M_INT(compton_xray_heating) M_DBL(reference_dust_to_gas_ratio)
  M_INT(cmb_temperature_floor) M_INT(agora_cmb_temperature_floor)
  M_DBL(redshift) M_INT(background_mode) M_INT(ionising_background_model)
  M_INT(lw_background_model) M_INT(isrf_background_model)
  M_DBL(isrf_background_habing) M_INT(dust_background_attenuation)
  M_STR(uv_rates_file) M_STR(ionising_rates_file) M_STR(lw_rates_file)
  M_STR(isrf_rates_file) M_DBL(lw_background_j21) M_INT(isrf_from_field)
  M_DBL(dust_cross_section_lw_cgs) M_DBL(dust_cross_section_pe_cgs)
  M_DBL(cosmic_ray_initial_value) M_INT(rt_use_fields) M_DBL(rt_heating_rate)
  M_DBL(rt_hi_ionisation_rate) M_DBL(rt_hei_ionisation_rate)
  M_DBL(rt_heii_ionisation_rate) M_DBL(rt_h2_dissociation_rate)
  M_DBL(volumetric_heating_rates) M_DBL(specific_heating_rates)
  M_DBL(HydrogenFractionByMass) M_DBL(initial_hii_to_h_ratio)
  M_DBL(initial_hei_to_h_ratio) M_DBL(initial_heii_to_h_ratio)
  M_DBL(initial_heiii_to_h_ratio) M_DBL(initial_di_to_h_ratio)
  M_DBL(initial_dii_to_h_ratio) M_DBL(initial_hm_to_h_ratio)
  M_DBL(initial_h2i_to_h_ratio) M_DBL(initial_h2ii_to_h_ratio)
  M_DBL(initial_hdi_to_h_ratio) M_INT(uv_self_shielding)
  M_DBL(uv_shielding_threshold) M_INT(h2_shielding_length_method)
  M_DBL(h2_shielding_path_in_kernel_radii) M_INT(uv_shielding_length_method)
  M_DBL(uv_shielding_length_factor) M_DBL(uv_shielding_path_in_kernel_radii)
  M_STR(uv_shielding_table_file) M_INT(max_substeps) M_INT(subcycle_on_failure)
  M_DBL(energy_change_limit) M_INT(cheap_energy_update) M_INT(solver)
  M_DBL(accurate_relative_tolerance) M_DBL(accurate_absolute_tolerance)
  M_DBL(thermal_time) M_INT(table_density_edge) M_DBL(maximal_density)
  printf("MOD cosmic_ray_context=%d\n", c.cosmic_ray_context != NULL);
  printf("KEYPASS_DONE\n");
  return 0;
}
