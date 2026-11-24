// One-zone thermal equilibrium at fixed density with the cooling library that
// SWIFT links, with the flags SWIFT sets for the ISRF: dust chemistry, a
// per-particle field G0 and no UV background.
//
// Build with Grackle:
//   gcc -O1 onezone.c -o onezone_grackle -I$GRACKLE_ROOT/include
//       -L$GRACKLE_ROOT/lib -lgrackle -lhdf5 -lm -Wl,-rpath,$GRACKLE_ROOT/lib
// Build with geajay:
//   gcc -O1 -DUSE_GEAJAY onezone.c -o onezone_geajay -I$GEAJAY_ROOT/include
//       -L$GEAJAY_ROOT/lib -lgeajay -lhdf5 -lm -Wl,-rpath,$GEAJAY_ROOT/lib
//
// Usage: onezone <network> <pe_model> <n_H [cm^-3]> <G0 [Habing]>
//                <Z [mass fraction]> <T0 [K]> <x_e0> <t_end [Myr]>
//                <cloudy table> [dt_max [Myr]]
// The last line starts with FINAL and gives the equilibrium temperature,
// specific energy and relaxation time.
#define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_GEAJAY
#include <geajay.h>
#include <geajay_float.h>
typedef geajay_code_units code_units_t;
typedef geajay_chemistry_data chem_t;
typedef geajay_chemistry_data_storage rates_t;
typedef geajay_field_data fields_t;
typedef geajay_float real_t;
#define set_defaults geajay_global_set_default_chemistry_parameters
#define init_chem geajay_initialize_chemistry_data
#define init_fields geajay_initialize_field_data
#define set_vel geajay_set_velocity_units
#define solve geajay_solve_chemistry
#define cooling_time geajay_calculate_cooling_time
#define temperature geajay_calculate_temperature
#define pressure geajay_calculate_pressure
#else
#include <grackle.h>
typedef code_units code_units_t;
typedef chemistry_data chem_t;
typedef chemistry_data_storage rates_t;
typedef grackle_field_data fields_t;
typedef gr_float real_t;
#define set_defaults set_default_chemistry_parameters
#define init_chem(c, r, u) local_initialize_chemistry_data(c, r, u)
#define init_fields(f) memset(f, 0, sizeof(*(f)))
#define set_vel set_velocity_units
#define solve(c, r, u, f, dt) local_solve_chemistry(c, r, u, f, dt)
#define cooling_time(c, r, u, f, t) local_calculate_cooling_time(c, r, u, f, t)
#define temperature(c, r, u, f, t) local_calculate_temperature(c, r, u, f, t)
#define pressure(c, r, u, f, p) local_calculate_pressure(c, r, u, f, p)
#endif

static const double mp = 1.67262192e-24, kb = 1.380649e-16, X = 0.76,
                    Myr = 3.15576e13;

int main(int argc, char **argv) {
  if (argc < 10) return 2;
  const int network = atoi(argv[1]), pe = atoi(argv[2]);
  const double nH = atof(argv[3]), G0 = atof(argv[4]), Z = atof(argv[5]),
               T0 = atof(argv[6]);
  const double xe0 = atof(argv[7]), t_end = atof(argv[8]);
  const double dt_max = argc > 10 ? atof(argv[10]) : 0.01 / nH;
  const double Y = 1. - X;

  code_units_t u;
  u.comoving_coordinates = 0;
  u.density_units = mp;
  u.length_units = 3.0856775814913673e21;
  u.time_units = Myr;
  u.a_units = 1.;
  u.a_value = 1.;
  set_vel(&u);
  const double vel2 = pow(u.length_units / u.time_units, 2);

  chem_t *c = malloc(sizeof(chem_t));
  if (!set_defaults(c)) return 3;
  c->use_grackle = 1;
  c->with_radiative_cooling = 1;
  c->primordial_chemistry = network;
  c->metal_cooling = 1;
  c->UVbackground = 0;
  c->dust_chemistry = 1;
  c->use_isrf_field = 1;
  c->photoelectric_heating = pe;
  c->HydrogenFractionByMass = X;
  c->cmb_temperature_floor = 1;
  c->CaseBRecombination = 1;
  c->max_iterations = 10000;
  c->exit_after_iterations_exceeded = 1;
  c->grackle_data_file = argv[9];
  rates_t r;
  if (!init_chem(c, &r, &u)) return 4;

  const double rho = nH * mp / X / u.density_units;
  real_t d = rho, e, vx = 0, vy = 0, vz = 0, isrf = G0, metal = Z * rho;
  real_t HI = X * rho * (1 - xe0), HII = X * rho * xe0, HeI = Y * rho, HeII = 0,
         HeIII = 0;
  real_t de = xe0 * X * rho, zero = 0;
  const double mu = 1. / (X * (1 + xe0) + Y / 4);
  e = kb * T0 / ((5. / 3 - 1) * mu * mp) / vel2;
  fields_t f;
  init_fields(&f);
  int dims[3] = {1, 1, 1}, ss[3] = {0, 0, 0}, ee[3] = {0, 0, 0};
  f.grid_rank = 3;
  f.grid_dimension = dims;
  f.grid_start = ss;
  f.grid_end = ee;
  f.density = &d;
  f.internal_energy = &e;
  f.x_velocity = &vx;
  f.y_velocity = &vy;
  f.z_velocity = &vz;
  f.metal_density = &metal;
  f.isrf_habing = &isrf;
  if (network >= 1) {
    f.HI_density = &HI;
    f.HII_density = &HII;
    f.HeI_density = &HeI;
    f.HeII_density = &HeII;
    f.HeIII_density = &HeIII;
    f.e_density = &de;
  }
  if (network >= 2) {
    f.HM_density = &zero;
    f.H2I_density = &zero;
    f.H2II_density = &zero;
  }
  if (network >= 3) {
    f.DI_density = &zero;
    f.DII_density = &zero;
    f.HDI_density = &zero;
  }

  real_t tc = 0, T = 0, P = 0;
  double t = 0;
  while (t < t_end) {
    if (!cooling_time(c, &r, &u, &f, &tc)) return 5;
    double dt = fmin(0.1 * fabs((double)tc), dt_max);
    if (!(dt > 0)) dt = dt_max;
    if (!solve(c, &r, &u, &f, dt)) return 6;
    t += dt;
  }
  temperature(c, &r, &u, &f, &T);
  pressure(c, &r, &u, &f, &P);
  const double T_end = T, u_end = e;

  /* Relaxation time at the end state with the species frozen: tau = -de /
   * d(de/dt). */
  cooling_time(c, &r, &u, &f, &tc);
  const double rate0 = e / (double)tc; /* de/dt: tc is positive for heating */
  const double e0 = e;
  e = e0 * 1.001;
  cooling_time(c, &r, &u, &f, &tc);
  const double rate1 = e / (double)tc;
  e = e0;
  const double tau =
      (rate0 - rate1) > 0 ? (0.001 * e0) / (rate0 - rate1) : INFINITY;

  /* Drift over a further tenth of the elapsed time. */
  solve(c, &r, &u, &f, 0.1 * t);
  temperature(c, &r, &u, &f, &T);
  printf(
      "FINAL n_H=%g G0=%g pe=%d T0=%g T=%.8g u_cgs=%.8g P_over_k=%.8g xe=%.4e "
      "tau_Myr=%.5g drift=%.3e\n",
      nH, G0, pe, T0, T_end, u_end * vel2,
      (double)P * u.density_units * vel2 / kb, (double)de / (X * rho), tau,
      ((double)T - T_end) / T_end);
  return 0;
}
