/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/**
 * @brief Particle fields not needed during the SPH loops over neighbours.
 *
 * This structure contains the particle fields that are not used in the
 * density or force loops. Quantities should be used in the kick, drift and
 * potentially ghost tasks only.
 */
struct xpart {

  float x_diff[3]; /*!< Offset between current position and position at last
                      tree rebuild. */

  float v_full[3]; /*!< Velocity at the last full step. */

  float u_full; /*!< Thermal energy at the last full step. */

} __attribute__((aligned(xpart_align)));

/**
 * @brief Particle fields for the SPH particles
 *
 * The density and force substructures are used to contain variables only used
 * within the density and force loops over neighbours. All more permanent
 * variables should be declared in the main part of the part structure,
 */
struct part {

  double x[3]; /*!< Particle position. */

  float v[3]; /*!< Particle predicted velocity. */

  float a_hydro[3]; /*!< Particle acceleration. */

  float mass; /*!< Particle mass. */

  float h; /*!< Particle smoothing length. */

  int ti_begin; /*!< Time at the beginning of time-step. */

  int ti_end; /*!< Time at the end of time-step. */

  float u; /*!< Particle internal energy. */

  float u_dt; /*!< Time derivative of the internal energy. */

  float rho; /*!< Particle density. */

  float rho_dh; /*!< Derivative of density with respect to h */

  /* Store density/force specific stuff. */
  union {

    /**
     * @brief Structure for the variables only used in the density loop over
     * neighbours.
     *
     * Quantities in this sub-structure should only be accessed in the density
     * loop over neighbours and the ghost task.
     */
    struct {

      float wcount; /*!< Neighbour number count. */

      float wcount_dh; /*!< Derivative of the neighbour number with respect to
                          h. */
    } density;

    /**
     * @brief Structure for the variables only used in the force loop over
     * neighbours.
     *
     * Quantities in this sub-structure should only be accessed in the force
     * loop over neighbours and the ghost, drift and kick tasks.
     */
    struct {

      float pressure; /*!< Particle pressure. */

      float v_sig; /*!< Particle signal velocity */

      float h_dt; /*!< Time derivative of smoothing length  */

    } force;
  };

  long long id; /*!< Particle unique ID. */

  struct gpart* gpart; /*!< Pointer to corresponding gravity part. */

} __attribute__((aligned(part_align)));
