#include "config.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

/* Local includes */
#include "active.h"
#include "adiabatic_index.h"
#include "dimension.h"
#include "engine.h"
#include "error.h"
#include "hydro.h"
#include "space.h"
#include "timeline.h"

/**
 * This file contains the routines for driven turbulence/stirring; use for
 * things like idealized turbulence tests, large-eddy simulations, and the like.
 *
 *  This code was originally written for GADGET3 by Andreas Bauer; it has been
 *   modified slightly by Phil Hopkins for GIZMO, but is largely intact.
 */

/***************************************/

/* Andreas Bauer's paper on turbulence:
 // sub-sonic (Mach~0.3) test: //
 ST_decay        1.
 ST_energy       0.0002 (sigma=0.014)
 ST_DtFreq       0.005
 ST_Kmin         6.27
 ST_Kmax         12.57
 ST_SolWeight    1.
 ST_AmplFac      1.
 ST_Seed         42
 ST_SpectForm    2

 // trans-sonic (Mach~1.2/3.5) test: //
 ST_decay        0.5
 ST_energy       0.21 (sigma=0.21-3.0)
 ST_DtFreq       0.005
 ST_Kmin         6.27
 ST_Kmax         12.57
 ST_SolWeight    1.
 ST_AmplFac      1.
 ST_Seed         42
 ST_SpectForm    2

 // super-sonic (Mach~8.4) test: //
 ST_decay        0.05
 ST_energy       25.0 (sigma=12.247)
 ST_DtFreq       0.005
 ST_Kmin         6.27
 ST_Kmax         18.85
 ST_SolWeight    1.
 ST_AmplFac      1.
 ST_Seed         42
 ST_SpectForm    1
 */

/***************************************/

const int TurbDriving_Global_DrivingSpectrumKey =
    0;  // driving pwr-spec: 0=Ek~const; 1=sharp-peak at kc; 2=Ek~k^(-5/3);
        // 3=Ek~k^-2
const int TurbDriving_Global_DrivingRandomNumberKey =
    1;  // random number seed for modes
const double TurbDriving_Global_DrivingScaleKMinVar =
    1.0;  // minimum driving-k: should be ~2.*M_PI/All.BoxSize -->
          // TurbDrive_MaxWavelength
const double TurbDriving_Global_DrivingScaleKMaxVar =
    2.0;  // maximum driving-k: set to couple times Kmin or more if more cascade
          // desired   // should be < MaxWavelength
const double TurbDriving_Global_AccelerationPowerVariable =
    1.0;  // energy of driving-scale modes: sets norm of turb -->
          // TurbDrive_ApproxRMSVturb
const double TurbDriving_Global_DecayTime =
    0.;  // decay time for driving-mode phase correlations   -->
         // TurbDrive_CoherenceTime
const double TurbDriving_Global_SolenoidalFraction =
    0.5;  // fractional wt of solenoidal modes (wt*curl + (1-wt)*div) -->
          // TurbDrive_SolenoidalFraction
const double TurbDriving_Global_DtTurbUpdates =
    1.0;  // time interval for driving updates (set by hand)  -->
          // TurbDrive_TimeBetweenTurbUpdates

/* block of global variables used specifically for the set of subroutines below,
 * which need to be carried between timesteps */
double *StOUPhases;     // random fluctuating component of the amplitudes
double *StAmpl;         // relative amplitude for each k
double *StAka;          // phases (real part)
double *StAkb;          // phases (imag part)
double *StMode;         // k vectors
int StNModes;           // total number of modes
integertime_t StTPrev;  // time of last update (to determine when next will be)
gsl_rng *StRng;         // random number generator key
FILE *FdTurb;
double TurbInjectedEnergy;
double TurbDissipatedEnergy;

/* routine to return gaussian random number with zero mean and unity variance */
double st_turbdrive_get_gaussian_random_variable(void) {
  double r0 = gsl_rng_uniform(StRng), r1 = gsl_rng_uniform(StRng);
  return sqrt(2. * log(1. / r0)) * cos(2. * M_PI * r1);
}

/* return the driving scale needed for scaling some other quantities below,
 * corresponding to our global variable convention */
double st_return_driving_scale(void) {
  return TurbDriving_Global_DrivingScaleKMinVar;  // this is now spatial scale
}

/* return the coherence time of the driving scale modes. if global variable
 * negative, it uses the eddy turnover time of the driving-scale modes */
double st_return_mode_correlation_time(void) {
  if (TurbDriving_Global_DecayTime > 0)
    return TurbDriving_Global_DecayTime;
  else
    return st_return_driving_scale() /
           TurbDriving_Global_AccelerationPowerVariable;
}

/* return the rms acceleration we expect, using either the 'dissipation rate' or
 * 'turbulent velocity' conventions for our variables */
double st_return_rms_acceleration(void) {
  return TurbDriving_Global_AccelerationPowerVariable /
         st_return_mode_correlation_time();  // new convention, hoping this is
                                             // more clear re: meaning of
                                             // variable
}

/* initialize phase variables */
void st_turbdrive_init_ouseq(void) {
  for (int i = 0; i < 6 * StNModes; i++) {
    StOUPhases[i] = st_turbdrive_get_gaussian_random_variable() *
                    st_return_rms_acceleration();
  }
}

/* return time interval between turbulent driving field updates based on global
 * variable. if negative, default to small interval of coherence time by
 * default. */
double st_return_dt_between_updates(void) {
  if (TurbDriving_Global_DtTurbUpdates > 0)
    return TurbDriving_Global_DtTurbUpdates;
  else
    return 0.01 * st_return_mode_correlation_time();
}

/* update the Markov random variable that is the dimensional multiplier for the
 * acceleration field, which has a correlation time specified */
void st_update_ouseq(void) {
  double damping =
      exp(-st_return_dt_between_updates() / st_return_mode_correlation_time());
  for (int i = 0; i < 6 * StNModes; i++) {
    StOUPhases[i] = StOUPhases[i] * damping +
                    st_return_rms_acceleration() *
                        sqrt(1. - damping * damping) *
                        st_turbdrive_get_gaussian_random_variable();
  }
}

/* routine to calculate the projected phases/acceleration field variables, using
 * the fourier-space solenoidal/compressible projection */
void st_turbdrive_calc_phases(void) {
  int i, j;
  for (i = 0; i < StNModes; i++) {
    double ka = 0., kb = 0., kk = 0.;
    int dim = hydro_dimension;
    for (j = 0; j < dim; j++) {
      kk += StMode[3 * i + j] * StMode[3 * i + j];
      ka += StMode[3 * i + j] * StOUPhases[6 * i + 2 * j + 1];
      kb += StMode[3 * i + j] * StOUPhases[6 * i + 2 * j + 0];
    }
    for (j = 0; j < dim; j++) {
      double diva = StMode[3 * i + j] * ka / kk;
      double divb = StMode[3 * i + j] * kb / kk;
      double curla = StOUPhases[6 * i + 2 * j + 0] - divb;
      double curlb = StOUPhases[6 * i + 2 * j + 1] - diva;

      StAka[3 * i + j] = TurbDriving_Global_SolenoidalFraction * curla +
                         (1. - TurbDriving_Global_SolenoidalFraction) * divb;
      StAkb[3 * i + j] = TurbDriving_Global_SolenoidalFraction * curlb +
                         (1. - TurbDriving_Global_SolenoidalFraction) * diva;
    }
  }
}

/* parent routine to initialize and update turbulent driving fields and to track
 * different variables used for analyzing power spectra of dissipation, etc. */
void set_turb_ampl(const struct engine *e, struct cell *c) {
  double delta = (e->ti_current - StTPrev) * e->time_base;
  double Dt_Update = st_return_dt_between_updates();

  const int count = c->hydro.count;

  if (delta >= Dt_Update) {
    if (delta > 0) {
      double e_diss_sum = 0, e_drive_sum = 0, glob_diss_sum = 0,
             glob_drive_sum = 0;
      // message(" ..updating fields tracked for following injected energy and
      // dissipation");

      for (int i = 0; i < count; i++) {

        struct part *p = &c->hydro.parts[i];

        if (part_is_inhibited(p, e)) continue;

        const float m = hydro_get_mass(p);

        if (m > 0.) {
          e_diss_sum += p->turbulence.energy_diss;
          p->turbulence.du_dt_diss = (p->turbulence.energy_diss / m) / delta;
          p->turbulence.energy_diss = 0;

          e_drive_sum += p->turbulence.energy_drive;
          p->turbulence.du_dt_drive = (p->turbulence.energy_drive / m) / delta;
          p->turbulence.energy_drive = 0;
        } else {
          p->turbulence.energy_diss = 0.f;
          p->turbulence.du_dt_diss = 0.f;
          p->turbulence.energy_drive = 0.f;
          p->turbulence.du_dt_drive = 0.f;
        }
      }

      // MATTHIEU: Replace by a local atomic reduction
      MPI_Allreduce(&e_diss_sum, &glob_diss_sum, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      MPI_Allreduce(&e_drive_sum, &glob_drive_sum, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
      TurbDissipatedEnergy += glob_diss_sum;
      TurbInjectedEnergy += glob_drive_sum;
    }

    // MATTHIEU: Need to do this only once!!!
    message(" ..updating fourier-space phase information");
    st_update_ouseq();
    message(
        " ..calculating coefficients and phases following desired projection");
    st_turbdrive_calc_phases();
    StTPrev = StTPrev + Dt_Update / e->time_base;
    message(" ..updated turbulent stirring field at time %f",
            StTPrev * e->time_base);
  }
}

/* routine to initialize the different modes and their relative amplitudes and
 * other global variables needed for the turbulence driving routines */
void init_turb(const struct space *s, const struct engine *e) {
  int ikx, iky, ikz;
  double kx, ky, kz, k, ampl;

  double kmin = TurbDriving_Global_DrivingScaleKMinVar;
  double kmax = TurbDriving_Global_DrivingScaleKMaxVar;
  kmin = 2. * M_PI / TurbDriving_Global_DrivingScaleKMinVar;
  kmax =
      2. * M_PI /
      TurbDriving_Global_DrivingScaleKMaxVar;  // convert these from spatial
                                               // lengths to wavenumbers in the
                                               // new convention we are using

  const double boxSize_X = s->dim[0];
  const double boxSize_Y = s->dim[1];
  const double boxSize_Z = s->dim[2];

  const int ikxmax = boxSize_X * kmax / 2. / M_PI;
  const int ikymax = boxSize_Y * kmax / 2. / M_PI;
  const int ikzmax = boxSize_Z * kmax / 2. / M_PI;

  StNModes = 0;
  for (ikx = 0; ikx <= ikxmax; ikx++) {
    kx = 2. * M_PI * ikx / boxSize_X;
    for (iky = 0; iky <= ikymax; iky++) {
      ky = 2. * M_PI * iky / boxSize_Y;
      for (ikz = 0; ikz <= ikzmax; ikz++) {
        kz = 2. * M_PI * ikz / boxSize_Z;
        k = sqrt(kx * kx + ky * ky + kz * kz);
        if (k >= kmin && k <= kmax) {
#ifdef HYDRO_DIMENSION_1D
          StNModes += 1;
#endif
#ifdef HYDRO_DIMENSION_2D
          StNModes += 2;
#endif
#ifdef HYDRO_DIMENSION_3D
          StNModes += 4;
#endif
        }
      }
    }
  }

  message(
      "Initializing turbulent driving: max integer mode number ikx/iky/ikz = "
      "%d %d %d",
      ikxmax, ikymax, ikzmax);
  StMode = (double *)swift_malloc("StModes", StNModes * 3 * sizeof(double));
  StAka = (double *)swift_malloc("StAka", StNModes * 3 * sizeof(double));
  StAkb = (double *)swift_malloc("StAkb", StNModes * 3 * sizeof(double));
  StAmpl = (double *)swift_malloc("StAmpl", StNModes * sizeof(double));
  StOUPhases =
      (double *)swift_malloc("StOUPhases", StNModes * 6 * sizeof(double));
  double kc = 0.5 * (kmin + kmax), amin = 0.,
         amplitude_integrated_allmodes = 0.;
  StNModes = 0;

  for (ikx = 0; ikx <= ikxmax; ikx++) {
    kx = 2. * M_PI * ikx / boxSize_X;

    for (iky = 0; iky <= ikymax; iky++) {
      ky = 2. * M_PI * iky / boxSize_Y;

      for (ikz = 0; ikz <= ikzmax; ikz++) {
        kz = 2. * M_PI * ikz / boxSize_Z;

        k = sqrt(kx * kx + ky * ky + kz * kz);
        if (k >= kmin && k <= kmax) {
          if (TurbDriving_Global_DrivingSpectrumKey == 0) {
            ampl = 1.;  // uniform amplitude for all
          } else if (TurbDriving_Global_DrivingSpectrumKey == 1) {
            ampl = 4.0 * (amin - 1.0) / ((kmax - kmin) * (kmax - kmin)) *
                       ((k - kc) * (k - kc)) +
                   1.0;  // spike at kc = k_driving
          } else if (TurbDriving_Global_DrivingSpectrumKey == 2) {
            // ampl = pow(k/kmin, (1.-NUMDIMS)- 5./3. ); // because this is
            // E[vector_k] for NUMDIMS, need extra NUMDIMS-1 power term here
            ampl = pow(
                k / kmin,
                1. / 3. - 0 * 0.5 * hydro_dimension);  // this should scale as
                                                       // the acceleration per
                                                       // eddy, ~v^2/L, crudely
          } else if (TurbDriving_Global_DrivingSpectrumKey == 3) {
            // ampl = pow(k/kmin, (1.-NUMDIMS)- 2. ); // because this is
            // E[vector_k] for NUMDIMS, need extra NUMDIMS-1 power term here
            ampl =
                pow(k / kmin,
                    0. - 0 * 0.5 * hydro_dimension);  // this should scale as
                                                      // the acceleration per
                                                      // eddy, ~v^2/L, crudely
          } else {
            error("unknown spectral form");
          }

          StAmpl[StNModes] = ampl;
          StMode[3 * StNModes + 0] = kx;
          StMode[3 * StNModes + 1] = ky;
          StMode[3 * StNModes + 2] = kz;
          message(
              "  Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, "
              "ampl=%f",
              StNModes, ikx, iky, ikz, StMode[3 * StNModes + 0],
              StMode[3 * StNModes + 1], StMode[3 * StNModes + 2],
              StAmpl[StNModes]);
          amplitude_integrated_allmodes += ampl * ampl;
          StNModes++;

#if defined(HYDRO_DIMENSION_2D) || defined(HYDRO_DIMENSION_3D)
          if (ikx > 0 ||
              iky > 0)  // if both of these are zero, only non-degenerate modes
                        // are the +/- z modes [ensured below]
          {
            StAmpl[StNModes] = ampl;
            if (iky != 0) {
              StMode[3 * StNModes + 0] = kx;
            } else {
              StMode[3 * StNModes + 0] = -kx;
            }
            StMode[3 * StNModes + 1] = -ky;
            StMode[3 * StNModes + 2] = kz;
            message(
                "  Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, "
                "ampl=%f",
                StNModes, ikx, -iky, ikz, StMode[3 * StNModes + 0],
                StMode[3 * StNModes + 1], StMode[3 * StNModes + 2],
                StAmpl[StNModes]);
            amplitude_integrated_allmodes += ampl * ampl;
            StNModes++;
          }

#if defined(HYDRO_DIMENSION_3D)
          if ((iky > 0 || ikz > 0) &&
              (ikx > 0 ||
               ikz >
                   0))  // if both of these are zero, only non-degenerate modes
                        // are the +/- x or +/- y modes [already ensured above]
          {
            StAmpl[StNModes] = ampl;
            if (ikz != 0) {
              StMode[3 * StNModes + 0] = kx;
            } else {
              StMode[3 * StNModes + 0] = -kx;
            }
            StMode[3 * StNModes + 1] = ky;
            StMode[3 * StNModes + 2] = -kz;
            message(
                "  Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, "
                "ampl=%f",
                StNModes, ikx, iky, -ikz, StMode[3 * StNModes + 0],
                StMode[3 * StNModes + 1], StMode[3 * StNModes + 2],
                StAmpl[StNModes]);
            amplitude_integrated_allmodes += ampl * ampl;
            StNModes++;
          }

          if ((ikx > 0 || iky > 0) && (ikx > 0 || ikz > 0) &&
              (iky > 0 || ikz > 0))  // if both of these are zero, only
                                     // non-degenerate modes are +/- z or +/- y
                                     // or +/- x modes already handled above
          {
            StAmpl[StNModes] = ampl;
            if (ikz == 0 || iky == 0) {
              StMode[3 * StNModes + 0] = -kx;
            } else {
              StMode[3 * StNModes + 0] = kx;
            }
            StMode[3 * StNModes + 1] = -ky;
            StMode[3 * StNModes + 2] = -kz;
            message(
                "  Mode: %d, ikx=%d, iky=%d, ikz=%d, kx=%f, ky=%f, kz=%f, "
                "ampl=%f",
                StNModes, ikx, -iky, -ikz, StMode[3 * StNModes + 0],
                StMode[3 * StNModes + 1], StMode[3 * StNModes + 2],
                StAmpl[StNModes]);
            amplitude_integrated_allmodes += ampl * ampl;
            StNModes++;
          }
#endif  // 3D
#endif  // 2D or 3D
        }
      }
    }
  }

  for (int i = 0; i < StNModes; i++) {
    StAmpl[i] *= sqrt(1. / amplitude_integrated_allmodes);
  }  // normalize total driving amplitude across all modes here

  StTPrev = -1;  // mark some arbitrarily old time as last update of turb
                 // driving fields
  StRng = gsl_rng_alloc(gsl_rng_ranlxd1);  // allocate seed variables
  gsl_rng_set(StRng,
              TurbDriving_Global_DrivingRandomNumberKey);  // initialize seed
  for (int j = 0; j < 100; j++) {
    double tmp;
    tmp = st_turbdrive_get_gaussian_random_variable();
    tmp++;
  }                            // cycle past initial seed
  st_turbdrive_init_ouseq();   // initialize variable for phases
  st_turbdrive_calc_phases();  // initialize phases
  // MATTHIEU: CAll this on all cells!
  set_turb_ampl(
      e, &e->s->cells_top[0]);  // set initial amplitudes and calculate initial
                                // quantities needed for dissipation measures
  StTPrev =
      e->ti_current;  // mark current time as last update of turb driving fields
}

/* return factor needed to renormalize below based on fraction of power
 * projected out in our solenoidal projection, to return the correct
 * normalization for accelerations. */
double solenoidal_frac_total_weight_renormalization(void) {
#if defined(HYDRO_DIMENSION_3D)
  return sqrt(3.0 / 3.0) * sqrt(3.0) * 1.0 /
         sqrt(1.0 - 2.0 * TurbDriving_Global_SolenoidalFraction +
              3.0 * TurbDriving_Global_SolenoidalFraction *
                  TurbDriving_Global_SolenoidalFraction);
#elif defined(HYDRO_DIMENSION_2D)
  return sqrt(3.0 / 2.0) * sqrt(3.0) * 1.0 /
         sqrt(1.0 - 2.0 * TurbDriving_Global_SolenoidalFraction +
              2.0 * TurbDriving_Global_SolenoidalFraction *
                  TurbDriving_Global_SolenoidalFraction);
#elif defined(HYDRO_DIMENSION_1D)
  return sqrt(3.0 / 1.0) * sqrt(3.0) * 1.0 /
         sqrt(1.0 - 2.0 * TurbDriving_Global_SolenoidalFraction +
              1.0 * TurbDriving_Global_SolenoidalFraction *
                  TurbDriving_Global_SolenoidalFraction);
#endif
}

/* routine to actually calculate the turbulent acceleration 'driving field'
 * force on every resolution element */
void add_turb_accel(const struct engine *e, struct cell *c) {
  set_turb_ampl(e, c);

  const int count = c->hydro.count;

  double acc[3], fac_sol = 2. * solenoidal_frac_total_weight_renormalization();
  for (int i = 0; i < count; ++i) {
    struct part *p = &c->hydro.parts[i];

    if (part_is_active(p, e)) {

      double fx = 0, fy = 0, fz = 0;
      for (int m = 0; m < StNModes; m++)  // calc force
      {
        double kxx = StMode[3 * m + 0] * p->x[0];
        double kyy = StMode[3 * m + 1] * p->x[1];
        double kzz = StMode[3 * m + 2] * p->x[2];
        double kdotx = kxx + kyy + kzz, ampl = StAmpl[m], realt = cos(kdotx),
               imagt = sin(kdotx);

        fx += ampl * (StAka[3 * m + 0] * realt - StAkb[3 * m + 0] * imagt);
        fy += ampl * (StAka[3 * m + 1] * realt - StAkb[3 * m + 1] * imagt);
        fz += ampl * (StAka[3 * m + 2] * realt - StAkb[3 * m + 2] * imagt);
      }
      fx *= fac_sol;
      fy *= fac_sol;
      fz *= fac_sol;

      const float mass = hydro_get_mass(p);

      if (mass > 0.) {
        acc[0] = fx;
        acc[1] = fy;
        acc[2] = fz;

        for (int j = 0; j < 3; j++) {
          p->turbulence.a[j] = acc[j];
        }
      } else {
        p->turbulence.a[0] = p->turbulence.a[1] = p->turbulence.a[2] = 0;
      }
    }
  }
  message("Finished turbulence driving (acceleration) computation");
}

/* routine to integrate the turbulent driving forces, specifically the
 * 'TurbAccel' variables that need to be drifted and kicked: note that we
 * actually do drifting and kicking in the normal routines, this is just to
 * integrate the dissipation rates etc used for our tracking */
void do_turb_driving_step_first_half(struct cell *c, const struct engine *e) {

  const int count = c->hydro.count;

  for (int i = 0; i < count; ++i) {
    struct part *p = &c->hydro.parts[i];

    if (part_is_inhibited(p, e)) continue;
    if (!part_is_starting(p, e)) continue;

    const integertime_t ti_step = get_integer_timestep(p->time_bin);
    const double dt_kick = (ti_step / 2) * e->time_base;

    double mass = hydro_get_mass(p);
    double ekin0 = 0.5 * mass *
                   (p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]);
    double vtmp[3];
    double dvel[3];
    for (int j = 0; j < 3; j++) {
      dvel[j] = p->turbulence.a[j] * dt_kick;
      vtmp[j] = p->v[j] + dvel[j];
    }
    double ekin1 = 0.5 * mass *
                   (vtmp[0] * vtmp[0] + vtmp[1] * vtmp[1] + vtmp[2] * vtmp[2]);
    p->turbulence.energy_drive += ekin1 - ekin0;
  }
}

/* routine to integrate the turbulent driving forces, specifically the
 * 'TurbAccel' variables that need to be drifted and kicked: note that we
 * actually do drifting and kicking in the normal routines, this is just to
 * integrate the dissipation rates etc used for our tracking */
void do_turb_driving_step_second_half(struct cell *c, const struct engine *e) {
  const int count = c->hydro.count;

  for (int i = 0; i < count; ++i) {
    struct part *p = &c->hydro.parts[i];

    if (part_is_inhibited(p, e)) continue;
    if (!part_is_active(p, e)) continue;

    const integertime_t ti_step = get_integer_timestep(p->time_bin);
    const double dt_kick = (ti_step / 2) * e->time_base;

    double mass = hydro_get_mass(p);
    double ekin0 = 0.5 * mass *
                   (p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]);
    double vtmp[3];
    double dvel[3];
    for (int j = 0; j < 3; j++) {
      dvel[j] = p->turbulence.a[j] * dt_kick;
      vtmp[j] = p->v[j] + dvel[j];
    }
    double ekin1 = 0.5 * mass *
                   (vtmp[0] * vtmp[0] + vtmp[1] * vtmp[1] + vtmp[2] * vtmp[2]);
    p->turbulence.energy_drive += ekin1 - ekin0;
  }
}

/* routine to record and optionally write to output files various statistics of
 * driven turbulence here (most relevant to idealized sims with a hard-coded
 * adiabatic EOS */
void log_turb_temp(struct engine *e) {
  struct part *parts = e->s->parts;
  const int nr_parts = e->s->nr_parts;

  double dudt_drive = 0, dudt_diss = 0, mass = 0, ekin = 0, ethermal = 0;

  for (int i = 0; i < nr_parts; i++) {

    struct part *p = &parts[i];

    if (part_is_inhibited(p, e)) continue;

    const float m = hydro_get_mass(p);

    dudt_drive += m * p->turbulence.du_dt_drive;
    dudt_diss += m * p->turbulence.du_dt_diss;
    ekin +=
        0.5 * m * (p->v[0] * p->v[0] + p->v[1] * p->v[1] + p->v[2] * p->v[2]);
    ethermal += m * p->u;
    mass += m;
  }
  double mach = sqrt(2. * ekin / (hydro_gamma * (hydro_gamma - 1) * ethermal));

  fprintf(FdTurb, "%g %g %g %g %g %g %g\n", e->time, mach,
          (ekin + ethermal) / mass, dudt_drive / mass, dudt_diss / mass,
          TurbInjectedEnergy / mass, TurbDissipatedEnergy / mass);
  fflush(FdTurb);
}
