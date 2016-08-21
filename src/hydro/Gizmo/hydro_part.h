/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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

/* Extra particle data not needed during the computation. */
struct xpart {

  /* Offset between current position and position at last tree rebuild. */
  float x_diff[3];

  /* Velocity at the last full step. */
  float v_full[3];

  /* Old density. */
  float omega;

} __attribute__((aligned(xpart_align)));

/* Data of a single particle. */
struct part {

  /* Particle position. */
  double x[3];

  /* Particle predicted velocity. */
  float v[3];

  /* Particle acceleration. */
  float a_hydro[3];

  /* Particle cutoff radius. */
  float h;

  /* Particle time of beginning of time-step. */
  int ti_begin;

  /* Particle time of end of time-step. */
  int ti_end;

  /* The primitive hydrodynamical variables. */
  struct {

    /* Fluid velocity. */
    float v[3];

    /* Density. */
    float rho;

    /* Pressure. */
    float P;

    /* Gradients of the primitive variables. */
    struct {

      /* Density gradients. */
      float rho[3];

      /* Fluid velocity gradients. */
      float v[3][3];

      /* Pressure gradients. */
      float P[3];

    } gradients;

    /* Quantities needed by the slope limiter. */
    struct {

      /* Extreme values of the density among the neighbours. */
      float rho[2];

      /* Extreme values of the fluid velocity among the neighbours. */
      float v[3][2];

      /* Extreme values of the pressure among the neighbours. */
      float P[2];

      /* Maximal distance to all neighbouring faces. */
      float maxr;

    } limiter;

  } primitives;

  /* The conserved hydrodynamical variables. */
  struct {

    /* Fluid momentum. */
    float momentum[3];

    /* Fluid mass (this field already exists outside of this struct as well). */
    float mass;

    /* Fluid thermal energy (not per unit mass!). */
    float energy;

    /* Fluxes. */
    struct {

      /* Mass flux. */
      float mass;

      /* Momentum flux. */
      float momentum[3];

      /* Energy flux. */
      float energy;

    } flux;

  } conserved;

  /* Geometrical quantities used for hydro. */
  struct {

    /* Volume of the particle. */
    float volume;

    /* Geometrical shear matrix used to calculate second order accurate
       gradients */
    float matrix_E[3][3];

  } geometry;

  /* Variables used for timestep calculation (currently not used). */
  struct {

    /* Maximum fluid velocity among all neighbours. */
    float vmax;

  } timestepvars;

  /* Quantities used during the volume (=density) loop. */
  struct {

    /* Particle velocity divergence. */
    float div_v;

    /* Derivative of particle number density. */
    float wcount_dh;

    /* Particle velocity curl. */
    float curl_v[3];

    /* Particle number density. */
    float wcount;

  } density;

  /* Quantities used during the force loop. */
  struct {

    /* Needed to drift the primitive variables. */
    float h_dt;

    /* Physical time step of the particle. */
    float dt;

    /* Actual velocity of the particle. */
    float v_full[3];

  } force;

  /* Specific stuff for the gravity-hydro coupling. */
  struct {

    /* Previous value of the gravitational acceleration. */
    float old_a[3];

    /* Previous value of the mass flux vector. */
    float old_mflux[3];

    /* Current value of the mass flux vector. */
    float mflux[3];

  } gravity;

  /* Particle mass (this field is also part of the conserved quantities...). */
  float mass;

  /* Particle ID. */
  long long id;

  /* Associated gravitas. */
  struct gpart *gpart;

  /* Variables needed for the code to compile (should be removed/replaced). */
  float rho;

  /* Old internal energy flux */
  float du_dt;

} __attribute__((aligned(part_align)));
