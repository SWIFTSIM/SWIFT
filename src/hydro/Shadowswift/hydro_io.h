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

#include "adiabatic_index.h"
#include "equation_of_state.h"
#include "hydro.h"
#include "hydro_gradients.h"
#include "hydro_slope_limiters.h"
#include "io_properties.h"
#include "riemann.h"

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
                                UNIT_CONV_DENSITY, parts, primitives.rho);
}

/**
 * @brief Get the internal energy of a particle
 *
 * @param e #engine.
 * @param p Particle.
 * @return Internal energy of the particle
 */
INLINE static void convert_u(const struct engine* e, const struct part* p,
                             const struct xpart* xp, float* ret) {
  ret[0] = hydro_get_internal_energy(p);
}

/**
 * @brief Get the entropic function of a particle
 *
 * @param e #engine.
 * @param p Particle.
 * @return Entropic function of the particle
 */
INLINE static void convert_A(const struct engine* e, const struct part* p,
                             const struct xpart* xp, float* ret) {
  ret[0] = hydro_get_entropy(p);
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
#ifdef SHADOWFAX_TOTAL_ENERGY
  return p->conserved.energy;
#else
  if (p->conserved.mass > 0.) {
    float momentum2;

    momentum2 = p->conserved.momentum[0] * p->conserved.momentum[0] +
                p->conserved.momentum[1] * p->conserved.momentum[1] +
                p->conserved.momentum[2] * p->conserved.momentum[2];

    ret[0] = p->conserved.energy + 0.5f * momentum2 / p->conserved.mass;
  } else {
    ret[0] = 0.;
  }
#endif
}

INLINE static void convert_part_pos(const struct engine* e,
                                    const struct part* p,
                                    const struct xpart* xp, double* ret) {

  if (e->s->periodic) {
    ret[0] = box_wrap(p->x[0], 0.0, e->s->dim[0]);
    ret[1] = box_wrap(p->x[1], 0.0, e->s->dim[1]);
    ret[2] = box_wrap(p->x[2], 0.0, e->s->dim[2]);
  } else {
    ret[0] = p->x[0];
    ret[1] = p->x[1];
    ret[2] = p->x[2];
  }
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

  *num_fields = 13;

  /* List what we want to write */
  list[0] = io_make_output_field_convert_part(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, 1.f, parts, xparts,
      convert_part_pos, "Co-moving positions of the particles");

  list[1] = io_make_output_field(
      "Velocities", FLOAT, 3, UNIT_CONV_SPEED, 0.f, parts, primitives.v,
      "Peculiar velocities of the stars. This is (a * dx/dt) where x is the "
      "co-moving positions of the particles");

  list[2] = io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts,
                                 conserved.mass, "Masses of the particles");

  list[3] = io_make_output_field(
      "SmoothingLengths", FLOAT, 1, UNIT_CONV_LENGTH, 1.f, parts, h,
      "Co-moving smoothing lengths (FWHM of the kernel) of the particles");

  list[4] = io_make_output_field_convert_part(
      "InternalEnergies", FLOAT, 1, UNIT_CONV_ENERGY_PER_UNIT_MASS,
      -3.f * hydro_gamma_minus_one, parts, xparts, convert_u,
      "Co-moving thermal energies per unit mass of the particles");

  list[5] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           parts, id, "Unique IDs of the particles");

  list[6] = io_make_output_field("Accelerations", FLOAT, 3,
                                 UNIT_CONV_ACCELERATION, 1.f, parts, a_hydro,
                                 "Accelerations of the particles(does not "
                                 "work in non-cosmological runs).");

  list[7] = io_make_output_field("Densities", FLOAT, 1, UNIT_CONV_DENSITY, -3.f,
                                 parts, primitives.rho,
                                 "Co-moving mass densities of the particles");

  list[8] =
      io_make_output_field("Volumes", FLOAT, 1, UNIT_CONV_VOLUME, -3.f, parts,
                           cell.volume, "Co-moving volumes of the particles");

  list[9] = io_make_output_field("GradDensities", FLOAT, 3, UNIT_CONV_DENSITY,
                                 1.f, parts, primitives.gradients.rho,
                                 "Gradient densities of the particles");

  list[10] = io_make_output_field_convert_part(
      "Entropies", FLOAT, 1, UNIT_CONV_ENTROPY, 1.f, parts, xparts, convert_A,
      "Co-moving entropies of the particles");

  list[11] = io_make_output_field("Pressures", FLOAT, 1, UNIT_CONV_PRESSURE,
                                  -3.f * hydro_gamma, parts, primitives.P,
                                  "Co-moving pressures of the particles");

  list[12] = io_make_output_field_convert_part(
      "TotalEnergies", FLOAT, 1, UNIT_CONV_ENERGY, -3.f * hydro_gamma_minus_one,
      parts, xparts, convert_Etot, "Total (co-moving) energy of the particles");
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
}

/**
 * @brief Are we writing entropy in the internal energy field ?
 *
 * @return 1 if entropy is in 'internal energy', 0 otherwise.
 */
INLINE static int writeEntropyFlag(void) { return 0; }
