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
 * Usage: testGeajayKeys <parameter file>, testGeajayKeys --defaults, or
 * testGeajayKeys <parameter file> --particle <density> <energy> [<cosmic-ray
 * energy density to set>]: also pass one gas particle (internal units) through
 * the module and print the cosmic-ray value that reaches the library
 * (testGeajayCosmicRays.sh). With -DNO_CR_MODES the program builds against a
 * module without the cosmic-ray source modes. */
#define _GNU_SOURCE

/* Config parameters. */
#include <config.h>

#include <dlfcn.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "swift.h"

#include "cooling_io.h"
#include "io_properties.h"

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
  printf("LIB cr.coulomb_loss_coefficient=%.10g\n",
         s->coulomb_loss_coefficient);
  printf("LIB cr.hadronic_loss_coefficient=%.10g\n",
         s->hadronic_loss_coefficient);
  return real(context, s, error);
}

/* The per-particle entry points: the value that the module passes. */
static void print_particle_call(const char *name, double value,
                                const geajay_particle_inputs *in) {
  printf("PARTCALL %s cosmic_ray_value=%.17g density=%.17g metal=%.17g\n", name,
         value, in->density, in->metal);
}

int geajay_particle_solve_cosmic_ray_reduced(
    const geajay_cosmic_ray_reduced_context *context,
    const geajay_species_layout *layout, void *state,
    const geajay_particle_inputs *inputs, double cosmic_ray_value, double dt,
    int n_substeps, geajay_cosmic_ray_output *output,
    geajay_particle_status *status) {
  static int (*real)(const geajay_cosmic_ray_reduced_context *,
                     const geajay_species_layout *, void *,
                     const geajay_particle_inputs *, double, double, int,
                     geajay_cosmic_ray_output *, geajay_particle_status *);
  if (!real)
    real = dlsym(RTLD_NEXT, "geajay_particle_solve_cosmic_ray_reduced");
  print_particle_call("solve", cosmic_ray_value, inputs);
  /* GEAJAY_KEYS_FAIL_SOLVE=1: report a solve that failed after all retries. */
  if (getenv("GEAJAY_KEYS_FAIL_SOLVE") != NULL) {
    status->reason = GEAJAY_PARTICLE_ITERATIONS_EXCEEDED;
    printf("LIBOUT solve failed\n");
    return GEAJAY_FAIL;
  }
  const int ret = real(context, layout, state, inputs, cosmic_ray_value, dt,
                       n_substeps, output, status);
  printf("LIBOUT cosmic_ray_energy_loss_rate=%.9g\n",
         output->cosmic_ray_energy_loss_rate);
  return ret;
}

int geajay_particle_cooling_time_cosmic_ray_reduced(
    const geajay_cosmic_ray_reduced_context *context,
    const geajay_species_layout *layout, const void *state,
    const geajay_particle_inputs *inputs, double cosmic_ray_value,
    double *cooling_time) {
  static int (*real)(const geajay_cosmic_ray_reduced_context *,
                     const geajay_species_layout *, const void *,
                     const geajay_particle_inputs *, double, double *);
  if (!real)
    real = dlsym(RTLD_NEXT, "geajay_particle_cooling_time_cosmic_ray_reduced");
  print_particle_call("cooling_time", cosmic_ray_value, inputs);
  return real(context, layout, state, inputs, cosmic_ray_value, cooling_time);
}

int geajay_particle_energy_rates_cosmic_ray_reduced(
    const geajay_cosmic_ray_reduced_context *context,
    const geajay_species_layout *layout, const void *state,
    const geajay_particle_inputs *inputs, double cosmic_ray_value,
    geajay_particle_rates *rates) {
  static int (*real)(const geajay_cosmic_ray_reduced_context *,
                     const geajay_species_layout *, const void *,
                     const geajay_particle_inputs *, double,
                     geajay_particle_rates *);
  if (!real)
    real = dlsym(RTLD_NEXT, "geajay_particle_energy_rates_cosmic_ray_reduced");
  print_particle_call("energy_rates", cosmic_ray_value, inputs);
  return real(context, layout, state, inputs, cosmic_ray_value, rates);
}

/* The cosmic-ray snapshot fields of the module: name, units and the value that
 * the snapshot would hold for the particle. */
static void print_fields(const char *when, struct unit_system *us,
                         struct phys_const *phys_const,
                         struct cooling_function_data *c, struct part *p,
                         struct xpart *xp) {
  struct engine e;
  struct io_props list[64];
  memset(&e, 0, sizeof(e));
  e.internal_units = us;
  e.physical_constants = phys_const;
  e.cooling_func = c;
  const int n = cooling_write_particles(p, xp, list);
  for (int i = 0; i < n; i++) {
    if (strncmp(list[i].name, "CosmicRay", 9) != 0 &&
        strcmp(list[i].name, "DustToGasRatios") != 0)
      continue;
    float value = -1.f;
    if (list[i].convert_part_f) list[i].convert_part_f(&e, p, xp, &value);
    printf("FIELD %s %s dimension=%d units=%d value=%.9g\n", when, list[i].name,
           list[i].dimension, (int)list[i].units, value);
  }
  printf("FIELD %s count=%d\n", when, n);
}

/* The species mass fractions that the module sets at the start. */
static void print_species(const char *when, const struct xpart *xp) {
  const struct cooling_xpart_data *d = &xp->cooling_data;
#define SP(n) printf("SPECIES %s %s=%.9g\n", when, #n, (double)d->n);
  SP(HI_frac) SP(HII_frac) SP(HeI_frac) SP(HeII_frac) SP(HeIII_frac)
  SP(e_frac)
#if COOLING_GEAJAY_MODE >= 2
  SP(HM_frac) SP(H2I_frac) SP(H2II_frac)
#endif
#if COOLING_GEAJAY_MODE >= 3
  SP(DI_frac) SP(DII_frac) SP(HDI_frac)
#endif
#if COOLING_GEAJAY_MODE >= 4
  SP(DM_frac) SP(HDII_frac) SP(HeHII_frac)
#endif
#undef SP
}

/* A gas particle at the given density and energy, solved through the module:
 * what the module passes to the library is printed by the entry points above.
 * With CR_SET the energy density is first written through the documented
 * function of the module. */
static void run_particle(struct swift_params *params, struct unit_system *us,
                         struct phys_const *phys_const,
                         struct hydro_props *hydro_properties,
                         struct cooling_function_data *c, double rho, double u,
                         int set_value, double value) {
  struct part p;
  struct xpart xp;
  struct cosmology cosmo;
  struct pressure_floor_props pressure_floor;
  struct entropy_floor_properties floor_props;
  memset(&p, 0, sizeof(p));
  memset(&xp, 0, sizeof(xp));
  cosmology_init_no_cosmo(&cosmo);
  /* GEAJAY_KEYS_SCALE_FACTOR: the particle density is then comoving and the
   * energy comoving, as in a cosmological run. */
  const char *scale_factor = getenv("GEAJAY_KEYS_SCALE_FACTOR");
  if (scale_factor != NULL) {
    const double a = atof(scale_factor);
    cosmo.a = a;
    cosmo.z = 1. / a - 1.;
    cosmo.a_inv = 1. / a;
    cosmo.a2_inv = 1. / (a * a);
    cosmo.a3_inv = 1. / (a * a * a);
    cosmo.a_factor_internal_energy = 1. / (a * a);
    /* GEAJAY_KEYS_MEAN_BARYON_DENSITY: mean baryonic density today [g cm^-3].
     * Without it the run has none, as without cosmology. */
    const char *rho_b = getenv("GEAJAY_KEYS_MEAN_BARYON_DENSITY");
    if (rho_b != NULL)
      cosmo.mean_density_Omega_b =
          atof(rho_b) * cosmo.a3_inv / c->units.density_units;
  }
  pressure_floor_init(&pressure_floor, phys_const, us, hydro_properties,
                      params);
  entropy_floor_init(&floor_props, phys_const, us, hydro_properties, params);
  cooling_update(phys_const, &cosmo, &pressure_floor, c, NULL, 0.);

  p.id = 1;
  p.mass = 1.;
  p.h = 0.1f;
  p.rho = (float)rho;
  p.u = (float)u;
  /* GEAJAY_KEYS_METAL_FRACTION: the total metal mass fraction of the
   * particle, 0.01 without it. */
  const char *metal_fraction = getenv("GEAJAY_KEYS_METAL_FRACTION");
  p.chemistry_data
      .smoothed_metal_mass_fraction[GEAR_CHEMISTRY_ELEMENT_COUNT - 1] =
      metal_fraction != NULL ? (float)atof(metal_fraction) : 0.01f;
  cooling_first_init_part(phys_const, us, hydro_properties, NULL, c, &p, &xp);
  printf("PART metal_mass_fraction=%.17g\n",
         (double)chemistry_get_total_metal_mass_fraction_for_cooling(&p));
  printf("PART stored_value_after_first_init=%.17g\n",
         (double)xp.cooling_data.cosmic_ray_value);
  print_species("after_first_init", &xp);
  print_fields("after_first_init", us, phys_const, c, &p, &xp);

#ifndef NO_CR_MODES
  if (set_value) {
    cooling_set_cosmic_ray_energy_density(&xp, value);
    printf("PART stored_value_after_set=%.17g\n",
           cooling_get_cosmic_ray_energy_density(&xp));
    print_fields("after_set", us, phys_const, c, &p, &xp);
  }
#endif

  const double time_step = 1.e-6;
  double rates[GEAJAY_RATE_NCHANNELS_COSMIC_RAY];
  cooling_time(phys_const, us, hydro_properties, &cosmo, c, &p, &xp);
  cooling_energy_rates(phys_const, us, hydro_properties, &cosmo, c, &p, &xp,
                       rates);
  cooling_cool_part(phys_const, us, &cosmo, hydro_properties, &floor_props,
                    &pressure_floor, c, &p, &xp, time_step, time_step, 0.);
#ifndef NO_CR_MODES
  printf("PART loss_rate_getter=%.9g\n",
         cooling_get_cosmic_ray_energy_loss_rate(&xp));
  printf("PART stored_value_after_solve=%.17g\n",
         cooling_get_cosmic_ray_energy_density(&xp));
  /* The depletion of the documented recipe over 1e10 s. */
  if (cooling_cosmic_ray_source_is_per_particle()) {
    const double loss = cooling_get_cosmic_ray_energy_loss_rate(&xp);
    const double e = cooling_get_cosmic_ray_energy_density(&xp);
    cooling_set_cosmic_ray_energy_density(&xp, fmax(0., e - loss * 1.e10));
    printf("PART stored_value_after_recipe=%.17g\n",
           cooling_get_cosmic_ray_energy_density(&xp));
  }
#endif
  printf("PART loss_rate_field=%.9g\n",
         (double)xp.cooling_data.cosmic_ray_energy_loss_rate);
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
  M_DBL(initial_hdi_to_h_ratio) M_DBL(initial_dm_to_h_ratio)
  M_DBL(initial_hdii_to_h_ratio) M_DBL(initial_hehii_to_h_ratio)
  M_INT(uv_self_shielding)
  M_DBL(uv_shielding_threshold) M_INT(h2_shielding_length_method)
  M_DBL(h2_shielding_path_in_kernel_radii) M_INT(uv_shielding_length_method)
  M_DBL(uv_shielding_length_factor) M_DBL(uv_shielding_path_in_kernel_radii)
  M_STR(uv_shielding_table_file) M_INT(max_substeps) M_INT(subcycle_on_failure)
  M_DBL(energy_change_limit) M_INT(with_radiative_cooling)
  M_INT(cheap_energy_update) M_INT(solver)
  M_DBL(accurate_relative_tolerance) M_DBL(accurate_absolute_tolerance)
  M_DBL(thermal_time) M_INT(table_density_edge) M_DBL(maximal_density)
  printf("MOD cosmic_ray_context=%d\n", c.cosmic_ray_context != NULL);
  printf("MOD proton_mass_cgs=%.17g\n", c.proton_mass_cgs);
  printf("MOD local_dust_to_gas_ratio=%.17g\n",
         c.chemistry_data.local_dust_to_gas_ratio);
  printf("MOD SolarMetalFractionByMass=%.17g\n",
         c.chemistry_data.SolarMetalFractionByMass);
#ifndef NO_CR_MODES
  M_INT(cosmic_ray_source)
  M_DBL(cosmic_ray_fire2_energy_density)
#endif
  printf("MOD attenuation_dust_ratio_at_z0.001=%.10g\n",
         cooling_attenuation_dust_to_gas_ratio(&c, 0.001));
  printf("MOD density_units=%.17g\n", c.units.density_units);
  if (argc >= 5 && strcmp(argv[2], "--particle") == 0)
    run_particle(params, &us, &phys_const, &hydro_properties, &c, atof(argv[3]),
                 atof(argv[4]), argc >= 6, argc >= 6 ? atof(argv[5]) : 0.);
  printf("KEYPASS_DONE\n");
  return 0;
}
