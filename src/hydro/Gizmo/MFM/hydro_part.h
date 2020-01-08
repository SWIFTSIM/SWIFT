/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2014 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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
#ifndef SWIFT_GIZMO_MFM_HYDRO_PART_H
#define SWIFT_GIZMO_MFM_HYDRO_PART_H

/* Data of a single particle. */
struct part {

  /* Particle ID. */
  long long id;

  /* Associated gravitas. */
  struct gpart *gpart;

  /* Particle position. */
  double x[3];

  /* In MFM, the particle and fluid velocities are the same.
     We use an anonymous union to make sure we can reference the
     same array with both names. */
  union {
    /* Particle predicted velocity. */
    float v[3];

    /* Fluid velocity. */
    float fluid_v[3];
  };

  /* Particle acceleration. */
  float a_hydro[3];

  /* Particle smoothing length. */
  float h;

  /* Density. */
  float rho;

  /* Pressure. */
  float P;

  /* Entropic function. */
  float A;

  union {
    /* Quantities used during the volume (=density) loop. */
    struct {

      /* Derivative of particle number density. */
      float wcount_dh;

      /* Particle number density. */
      float wcount;

    } density;

    /* Quantities needed by the slope limiter. */
    struct {

      /* Extreme values of the density among the neighbours. */
      float rho[2];

      /* Extreme values of the fluid velocity among the neighbours. */
      float v[3][2];

      /* Extreme values of the pressure among the neighbours. */
      float P[2];

      /* Extreme values of the entropic function among the neighbours. */
      float A[2];

      /* Maximal distance to all neighbouring faces. */
      float maxr;

    } limiter;

    struct {
      /* Fluxes. */
      struct {

        /* No mass flux, since it is always zero. */

        /* Momentum flux. */
        float momentum[3];

        /* Energy flux. */
        float energy;

      } flux;

      /* Variables used for timestep calculation. */
      struct {

        /* Maximum signal velocity among all the neighbours of the particle. The
         * signal velocity encodes information about the relative fluid
         * velocities
         * AND particle velocities of the neighbour and this particle, as well
         * as
         * the sound speed of both particles. */
        float vmax;

      } timestepvars;

      /* Quantities used during the force loop. */
      struct {

        /* Needed to drift the primitive variables. */
        float h_dt;

        /* Maximum kinetic energy of all neighbouring particles (in the rest
         * frame of this particle). Used for the entropy switch. */
        float Ekinmax;

      } force;
    };
  };

  /* Gradients of the primitive variables. */
  struct {

    /* Density gradients. */
    float rho[3];

    /* Fluid velocity gradients. */
    float v[3][3];

    /* Pressure gradients. */
    float P[3];

    /* Entropic function gradients. */
    float A[3];

  } gradients;

  /* The conserved hydrodynamical variables. */
  struct {

    /* Fluid mass */
    float mass;

    /* Fluid momentum. */
    float momentum[3];

    /* Fluid thermal energy (not per unit mass!). */
    float energy;

    /* Fluid entropy. */
    float entropy;

  } conserved;

  /* Geometrical quantities used for hydro. */
  struct {

    /* Volume of the particle. */
    float volume;

    /* Geometrical shear matrix used to calculate second order accurate
       gradients */
    float matrix_E[3][3];

    /* Correction factor for wcount. */
    float wcorr;

  } geometry;

  /* Chemistry information */
  struct chemistry_part_data chemistry_data;

  /*! Black holes information (e.g. swallowing ID) */
  struct black_holes_part_data black_holes_data;

  /* Time-step length */
  timebin_t time_bin;

  /*! Time-step limiter information */
  struct timestep_limiter_data limiter_data;

#ifdef GIZMO_FLAG_VARIABLE
  /* Flag variable containing diagnostic information about the particle. */
  int flag_variable;
#endif

#ifdef SWIFT_DEBUG_CHECKS

  /* Time of the last drift */
  integertime_t ti_drift;

  /* Time of the last kick */
  integertime_t ti_kick;

#endif

} SWIFT_STRUCT_ALIGN;

#endif /* SWIFT_GIZMO_MFM_HYDRO_PART_H */
