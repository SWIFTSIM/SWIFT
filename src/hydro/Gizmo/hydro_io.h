/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_HYDRO_IO_H
#define SWIFT_GIZMO_HYDRO_IO_H

#include "adiabatic_index.h"
#include "hydro.h"
#include "hydro_gradients.h"
#include "hydro_slope_limiters.h"
#include "io_properties.h"
#include "riemann.h"

/* Set the description of the particle movement. */
#if defined(GIZMO_FIX_PARTICLES)
#define GIZMO_PARTICLE_MOVEMENT "Fixed particles."
#else
#define GIZMO_PARTICLE_MOVEMENT "Particles move with flow velocity."
#endif

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void hydro_read_particles(struct part* parts,
                                        struct io_props* list,
                                        int* num_fields) {

  *num_fields = 8;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, parts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                parts, conserved.mass);
  list[3] = io_make_input_field("SmoothingLength", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, h);
  list[4] = io_make_input_field("InternalEnergy", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts,
                                conserved.energy);
  list[5] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, id);
  list[6] = io_make_input_field("Accelerations", FLOAT, 3, OPTIONAL,
                                UNIT_CONV_ACCELERATION, parts, a_hydro);
  list[7] = io_make_input_field("Density", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_DENSITY, parts, rho);
}

/**
 * @brief Get the internal energy of a particle
 *
 * @param e #engine.
 * @param p Particle.
 * @param ret (return) Internal energy of the particle
 */
INLINE static void convert_u(const struct engine* e, const struct part* p,
                             const struct xpart* xp, float* ret) {

  ret[0] = hydro_get_comoving_internal_energy(p, xp);
}

/**
 * @brief Get the entropic function of a particle
 *
 * @param e #engine.
 * @param p Particle.
 * @param ret (return) Entropic function of the particle
 */
INLINE static void convert_A(const struct engine* e, const struct part* p,
                             const struct xpart* xp, float* ret) {
  ret[0] = hydro_get_comoving_entropy(p, xp);
}

/**
 * @brief Get the total energy of a particle
 *
 * @param e #engine.
 * @param p Particle.
 * @return Total energy of the particle
 */
INLINE static void convert_Etot(const struct engine* e, const struct part* p,
                                const struct xpart* xp, float* ret) {
#ifdef GIZMO_TOTAL_ENERGY
  ret[0] = p->conserved.energy;
#else
  if (p->conserved.mass > 0.0f) {
    const float momentum2 =
        p->conserved.momentum[0] * p->conserved.momentum[0] +
        p->conserved.momentum[1] * p->conserved.momentum[1] +
        p->conserved.momentum[2] * p->conserved.momentum[2];

    ret[0] = p->conserved.energy + 0.5f * momentum2 / p->conserved.mass;
  } else {
    ret[0] = 0.0f;
  }
#endif
}

INLINE static void convert_part_pos(const struct engine* e,
                                    const struct part* p,
                                    const struct xpart* xp, double* ret) {
  const struct space* s = e->s;
  if (s->periodic) {
    ret[0] = box_wrap(p->x[0], 0.0, s->dim[0]);
    ret[1] = box_wrap(p->x[1], 0.0, s->dim[1]);
    ret[2] = box_wrap(p->x[2], 0.0, s->dim[2]);
  } else {
    ret[0] = p->x[0];
    ret[1] = p->x[1];
    ret[2] = p->x[2];
  }
  if (e->snapshot_use_delta_from_edge) {
    ret[0] = min(ret[0], s->dim[0] - e->snapshot_delta_from_edge);
    ret[1] = min(ret[1], s->dim[1] - e->snapshot_delta_from_edge);
    ret[2] = min(ret[2], s->dim[2] - e->snapshot_delta_from_edge);
  }
}

INLINE static void convert_part_vel(const struct engine* e,
                                    const struct part* p,
                                    const struct xpart* xp, float* ret) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cosmology* cosmo = e->cosmology;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  const integertime_t ti_beg = get_integer_time_begin(ti_current, p->time_bin);
  const integertime_t ti_end = get_integer_time_end(ti_current, p->time_bin);

  /* Get time-step since the last kick */
  float dt_kick_grav, dt_kick_hydro;
  if (with_cosmology) {
    dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_beg, ti_current);
    dt_kick_grav -=
        cosmology_get_grav_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
    dt_kick_hydro = cosmology_get_hydro_kick_factor(cosmo, ti_beg, ti_current);
    dt_kick_hydro -=
        cosmology_get_hydro_kick_factor(cosmo, ti_beg, (ti_beg + ti_end) / 2);
  } else {
    dt_kick_grav = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
    dt_kick_hydro = (ti_current - ((ti_beg + ti_end) / 2)) * time_base;
  }

  /* Extrapolate the velocites to the current time */
  hydro_get_drifted_velocities(p, xp, dt_kick_hydro, dt_kick_grav, ret);

  /* Conversion from internal units to peculiar velocities */
  ret[0] *= cosmo->a_inv;
  ret[1] *= cosmo->a_inv;
  ret[2] *= cosmo->a_inv;
}

INLINE static void convert_part_potential(const struct engine* e,
                                          const struct part* p,
                                          const struct xpart* xp, float* ret) {

  if (p->gpart != NULL)
    ret[0] = gravity_get_comoving_potential(p->gpart);
  else
    ret[0] = 0.f;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void hydro_write_particles(const struct part* parts,
                                         const struct xpart* xparts,
                                         struct io_props* list,
                                         int* num_fields) {

  *num_fields = 11;

  /* List what we want to write */
  list[0] = io_make_output_field_convert_part(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, parts, xparts,
      convert_part_pos, "Co-moving positions of the particles");

  list[1] = io_make_output_field_convert_part(
      "Velocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, parts, xparts,
      convert_part_vel,
      "Peculiar velocities of the stars. This is (a * dx/dt) where x is the "
      "co-moving positions of the particles");

  list[2] =
      io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts,
                           conserved.mass, "Co-moving masses of the particles");

  list[3] = io_make_output_field(
      "SmoothingLengths", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, parts, h,
      "Co-moving smoothing lengths (FWHM of the kernel) of the particles");

  list[4] = io_make_output_field_convert_part(
      "InternalEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      3.f * hydro_gamma_minus_one, parts, xparts, convert_u,
      "Co-moving thermal energies per unit mass of the particles");

  list[5] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           parts, id, "Unique IDs of the particles");

  list[6] = io_make_output_field("Densities", FLOAT, 1, UNIT_CONV_DENSITY, -3.f,
                                 parts, rho,
                                 "Co-moving mass densities of the particles");

  list[7] = io_make_output_field_convert_part(
      "Entropies", FLOAT, 1, UNIT_CONV_ENTROPY, 0.f, parts, xparts, convert_A,
      "Co-moving entropies of the particles");

  list[8] = io_make_output_field("Pressures", FLOAT, 1, UNIT_CONV_PRESSURE,
                                 -3.f * hydro_gamma, parts, P,
                                 "Co-moving pressures of the particles");

  list[9] = io_make_output_field_convert_part(
      "TotalEnergies", FLOAT, 1, UNIT_CONV_ENERGY, -3.f * hydro_gamma_minus_one,
      parts, xparts, convert_Etot, "Total (co-moving) energy of the particles");

  list[10] = io_make_output_field_convert_part(
      "Potentials", FLOAT, 1, UNIT_CONV_POTENTIAL, -1.f, parts, xparts,
      convert_part_potential, "Gravitational potentials of the particles");
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
INLINE static void hydro_write_flavour(hid_t h_grpsph) {
  /* Gradient information */
  io_write_attribute_s(h_grpsph, "Gradient reconstruction model",
                       HYDRO_GRADIENT_IMPLEMENTATION);

  /* Slope limiter information */
  io_write_attribute_s(h_grpsph, "Cell wide slope limiter model",
                       HYDRO_SLOPE_LIMITER_CELL_IMPLEMENTATION);
  io_write_attribute_s(h_grpsph, "Piecewise slope limiter model",
                       HYDRO_SLOPE_LIMITER_FACE_IMPLEMENTATION);

  /* Riemann solver information */
  io_write_attribute_s(h_grpsph, "Riemann solver type",
                       RIEMANN_SOLVER_IMPLEMENTATION);

  /* Particle movement information */
  io_write_attribute_s(h_grpsph, "Particle movement", GIZMO_PARTICLE_MOVEMENT);
}

/**
 * @brief Are we writing entropy in the internal energy field ?
 *
 * @return 1 if entropy is in 'internal energy', 0 otherwise.
 */
INLINE static int writeEntropyFlag(void) { return 0; }

#endif /* SWIFT_GIZMO_HYDRO_IO_H */
