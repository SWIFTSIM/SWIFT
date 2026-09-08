/* Standalone Grackle thermal-equilibrium solver harness. Companion to
 * verify_isrf_pe_equilibrium.py, which compiles and runs this file; not
 * part of the SWIFT build.
 *
 * Unlike verify_photoelectric_heating_rate_grackle_harness.c (which
 * isolates the photoelectric heating RATE at one fixed T, with
 * dust_chemistry hardcoded to 0), this harness answers a different
 * question: given a fixed (G0, n_H, Z', primordial_chemistry) tuple, what
 * temperature does Grackle's OWN net-heating-equals-net-cooling balance
 * settle at, evolving under exactly the Grackle configuration a real
 * with_photoelectric_heating=1 SWIFT run uses (dust_chemistry=1,
 * photoelectric_heating=2, use_isrf_field=1, real Cloudy metal-line
 * cooling) -- i.e. this harness DOES exercise the dust-recombination-
 * cooling code path the 2026-09-08 NaN fix touched, which the rate-only
 * harness provably never did.
 *
 * Units: same density_units=m_H trick as the sibling harness (dom=1
 * exactly; see that file's header for the derivation), so a code-unit
 * mass density numerically equals a cm^-3 number density for pure
 * hydrogen gas. HydrogenFractionByMass is left at the caller's real
 * value here (unlike the sibling harness's HydrogenFractionByMass=1.0
 * shortcut), since Grackle's own internal nH bookkeeping (rate_functions.c,
 * cool1d_multi_g.F) uses it directly for the tabulated-mode (primordial_
 * chemistry=0) heating/cooling rates, and this harness's whole point is to
 * match a real run's configuration, not simplify it away.
 *
 * Method (time-integration to equilibrium, not a T-bisection): starting
 * from a caller-supplied T_initial_K and a near-neutral primordial
 * species mix, repeatedly call calculate_cooling_time() to get the
 * current net-rate timescale, take a step dt = SAFETY_FACTOR*|cooling
 * time| (capped in growth relative to the previous step, since the
 * cooling time itself diverges on approach to equilibrium and an
 * unclamped step would overshoot), then solve_chemistry(dt) to evolve
 * BOTH the internal energy and the chemical species self-consistently.
 * This mirrors the standard Grackle usage pattern (see e.g. pygrackle's
 * evolve_constant_density), adapted to use an adaptively recomputed dt.
 * In practice the adaptive dt overshoots near equilibrium and settles
 * into a persistent few-percent limit cycle rather than damping to a
 * point (confirmed by tracing, HARNESS_TRACE=1 env var), so convergence
 * is judged on the time-weighted mean temperature of consecutive
 * BURN_IN_ITERATIONS-separated windows agreeing to CONVERGENCE_TOL, not
 * on the instantaneous temperature settling.
 *
 * Prints one "DATA,<T_initial_K>,<T_equilibrium_K>,<n_iterations>,
 * <t_elapsed_s>,<converged 0/1>" CSV line per invocation (one tuple per
 * process, unlike the sibling harness's internal Z' sweep loop -- the
 * Python driver invokes this binary once per extracted SWIFT particle).
 */

#include <grackle.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MH_CGS 1.67262171e-24 /* g, matches Grackle's own internal mh */

#define SAFETY_FACTOR 0.05
#define MAX_DT_GROWTH 1.2
#define DT_CEILING_S 1.0e16 /* ~300 Myr; far beyond any physical timescale
                               this problem can produce, just a guard
                               against a pathological single-step jump */
/* An adaptive dt sized off the instantaneous cooling time overshoots near
 * equilibrium (the cooling time itself diverges there) and settles into a
 * persistent few-percent limit-cycle oscillation rather than damping to a
 * point -- confirmed by tracing (HARNESS_TRACE=1). A step-to-step or
 * trailing-band tolerance therefore never triggers on an already-settled
 * state. Convergence is judged instead by comparing the time-weighted
 * mean temperature (the physically meaningful "average heating balances
 * average cooling" equilibrium value) of two consecutive windows: once
 * consecutive window means agree, the oscillation itself has stopped
 * drifting, which is what "dT/dt -> 0" means for a system with this kind
 * of numerical wobble on top of the true relaxation. */
#define BURN_IN_ITERATIONS 3000
#define STATS_WINDOW_ITERATIONS 3000
#define CONVERGENCE_TOL 5.0e-3
#define MAX_ITERATIONS 60000

int main(int argc, char **argv) {
  if (argc != 9 && argc != 10) {
    fprintf(stderr,
            "usage: %s <primordial_chemistry> <metal_cooling 0/1> "
            "<cloudy_table_path> <hydrogen_fraction_by_mass> <G0_habing> "
            "<n_H_cgs> <Zprime> <T_initial_K> [u_measured_cgs]\n",
            argv[0]);
    return 1;
  }

  const int primordial_chemistry = atoi(argv[1]);
  const int metal_cooling = atoi(argv[2]);
  const char *cloudy_table = argv[3];
  const double hydrogen_fraction_by_mass = atof(argv[4]);
  const double G0 = atof(argv[5]);
  const double n_H_cgs = atof(argv[6]);
  const double Zprime = atof(argv[7]);
  const double T_initial_K = atof(argv[8]);
  const int have_u_measured = (argc == 10);
  const double u_measured_cgs = have_u_measured ? atof(argv[9]) : 0.0;

  grackle_verbose = 0;

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
    return 1;
  }

  grackle_data->use_grackle = 1;
  grackle_data->with_radiative_cooling = 1;
  grackle_data->primordial_chemistry = primordial_chemistry;
  grackle_data->metal_cooling = metal_cooling;
  grackle_data->UVbackground = 0;
  grackle_data->grackle_data_file = strdup(cloudy_table);

  /* This is the real GEARFeedback:with_photoelectric_heating=1 branch
   * cooling.c's cooling_init_grackle() sets (dust_chemistry bundles
   * photoelectric heating, dust-recombination cooling, and H2-on-dust
   * formation under one switch) -- deliberately NOT simplified to
   * dust_chemistry=0 the way the sibling rate-only harness is, since
   * exercising the dust-recombination-cooling path is this harness's
   * whole point. */
  grackle_data->dust_chemistry = 1;
  grackle_data->photoelectric_heating = 2;
  grackle_data->use_isrf_field = 1;
  grackle_data->use_dust_density_field = 0; /* dust2gas = fgr * metallicity */
  grackle_data->h2_on_dust = 0;             /* GrackleCooling:H2_on_dust default */
  grackle_data->three_body_rate = 0;        /* GrackleCooling:H2_three_body_rate default */
  grackle_data->cie_cooling = 0;            /* GrackleCooling:H2_cie_cooling default */
  grackle_data->H2_self_shielding = 0;      /* GrackleCooling:H2_self_shielding default */
  grackle_data->self_shielding_method = 0;  /* GrackleCooling:self_shielding_method<=0 */
  grackle_data->cmb_temperature_floor = 1;  /* matches this example's params.yml */
  grackle_data->CaseBRecombination = 1;     /* cooling_init_grackle always sets this */
  grackle_data->HydrogenFractionByMass = hydrogen_fraction_by_mass;
  grackle_data->Gamma = 5.0 / 3.0;

  if (!initialize_chemistry_data(&my_units)) {
    fprintf(stderr, "initialize_chemistry_data failed\n");
    return 1;
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

  /* Near-neutral primordial-gas initial guess; the species network (at
   * primordial_chemistry>=1) and the internal energy both then relax
   * self-consistently to their joint equilibrium under repeated
   * solve_chemistry() calls below, so this starting guess does not need
   * to already be close to the true equilibrium ionization state. */
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

    /* metallicity(i) = metal/density/SolarMetalFractionByMass = Zprime,
     * matching cooling.c's chemistry_get_total_metal_mass_fraction_for_
     * cooling()*rho -> Grackle's own metal_density field. */
    my_fields.metal_density[i] =
        Zprime * grackle_data->SolarMetalFractionByMass * my_fields.density[i];

    my_fields.x_velocity[i] = 0.0;
    my_fields.y_velocity[i] = 0.0;
    my_fields.z_velocity[i] = 0.0;
    my_fields.internal_energy[i] = T_initial_K / temperature_units;
    my_fields.volumetric_heating_rate[i] = 0.0;
    my_fields.specific_heating_rate[i] = 0.0;
    my_fields.RT_HI_ionization_rate[i] = 0.0;
    my_fields.RT_HeI_ionization_rate[i] = 0.0;
    my_fields.RT_HeII_ionization_rate[i] = 0.0;
    my_fields.RT_H2_dissociation_rate[i] = 0.0;
    my_fields.RT_heating_rate[i] = 0.0;
    my_fields.isrf_habing[i] = G0;
  }

  gr_float *temperature = malloc(field_size * sizeof(gr_float));

  if (have_u_measured) {
    /* Report the temperature Grackle's OWN internal mu convention assigns
     * to a caller-supplied specific internal energy, at this same (n_H,
     * Z', mode) state -- this harness's density_units=m_H/length_units=1cm/
     * time_units=1s choice makes the code-unit specific energy numerically
     * equal to erg/g, so no conversion is needed before assigning it. Only
     * meaningful at primordial_chemistry=0: that tabulated mode's mu(T,Z)
     * comes from the Cloudy table, not from the species fields, so this
     * near-neutral placeholder species IC does not bias it; at
     * primordial_chemistry>=1, mu genuinely depends on the actual
     * ionization state, which this placeholder IC does not represent, so
     * the caller must not rely on this value there. */
    my_fields.internal_energy[0] = u_measured_cgs;
    if (!calculate_temperature(&my_units, &my_fields, temperature)) {
      fprintf(stderr, "calculate_temperature failed (measurement)\n");
      return 1;
    }
    printf("MEASURE,%.10e\n", (double)temperature[0]);
  }

  /* Calibrate internal_energy so the achieved starting T really is
   * T_initial_K (same one-pass linear correction as the sibling harness;
   * mu is ~constant with internal_energy at fixed composition). */
  my_fields.internal_energy[0] = T_initial_K / temperature_units;
  if (!calculate_temperature(&my_units, &my_fields, temperature)) {
    fprintf(stderr, "calculate_temperature failed (calibration)\n");
    return 1;
  }
  my_fields.internal_energy[0] *= T_initial_K / temperature[0];

  gr_float *cooling_time = malloc(field_size * sizeof(gr_float));

  double dt = 0.0;
  double t_elapsed = 0.0;
  double window_dt_sum = 0.0, window_T_dt_sum = 0.0;
  double previous_window_mean = -1.0;
  double T_equilibrium = 0.0;
  int converged = 0;
  int iter;
  const int trace = getenv("HARNESS_TRACE") != NULL;

  for (iter = 0; iter < MAX_ITERATIONS; iter++) {
    if (!calculate_cooling_time(&my_units, &my_fields, cooling_time)) {
      fprintf(stderr, "calculate_cooling_time failed at iteration %d\n", iter);
      return 1;
    }

    const double dt_target = SAFETY_FACTOR * fabs((double)cooling_time[0]);
    double dt_next = (iter == 0) ? dt_target : fmin(dt_target, dt * MAX_DT_GROWTH);
    if (dt_next > DT_CEILING_S) dt_next = DT_CEILING_S;
    if (!(dt_next > 0.0)) {
      /* cooling_time is exactly 0 (should not happen away from a
       * pathological state); fall back to a tiny nonzero step rather
       * than stalling forever on dt=0. */
      dt_next = 1.0e-30;
    }
    dt = dt_next;

    if (!solve_chemistry(&my_units, &my_fields, dt)) {
      fprintf(stderr, "solve_chemistry failed at iteration %d (dt=%.6e)\n",
              iter, dt);
      return 1;
    }
    t_elapsed += dt;

    if (!calculate_temperature(&my_units, &my_fields, temperature)) {
      fprintf(stderr, "calculate_temperature failed at iteration %d\n", iter);
      return 1;
    }
    const double T_now = (double)temperature[0];

    if (trace && (iter % 5000 == 0)) {
      fprintf(stderr, "trace iter=%d t=%.6e s dt=%.6e s T=%.6e K\n", iter,
              t_elapsed, dt, T_now);
    }

    if (iter >= BURN_IN_ITERATIONS) {
      window_dt_sum += dt;
      window_T_dt_sum += T_now * dt;
      const int iter_in_window = iter - BURN_IN_ITERATIONS + 1;
      if (iter_in_window % STATS_WINDOW_ITERATIONS == 0) {
        const double window_mean = window_T_dt_sum / window_dt_sum;
        if (trace) {
          fprintf(stderr, "  window mean T=%.6e K (dt_sum=%.6e s)\n",
                  window_mean, window_dt_sum);
        }
        if (previous_window_mean > 0.0 &&
            fabs(window_mean - previous_window_mean) / previous_window_mean <
                CONVERGENCE_TOL) {
          converged = 1;
          T_equilibrium = window_mean;
          break;
        }
        previous_window_mean = window_mean;
        T_equilibrium = window_mean; /* best available estimate so far */
        window_dt_sum = 0.0;
        window_T_dt_sum = 0.0;
      }
    }
  }

  printf("DATA,%.10e,%.10e,%d,%.10e,%d\n", T_initial_K, T_equilibrium,
        iter + 1, t_elapsed, converged);

  if (!converged) {
    fprintf(stderr,
            "WARNING: did not converge within %d iterations (last window "
            "mean T=%.6e K) -- reported T_equilibrium is NOT trustworthy\n",
            MAX_ITERATIONS, T_equilibrium);
  }

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
  free((char *)my_chemistry->grackle_data_file);
  free(my_chemistry);

  return converged ? 0 : 2;
}
