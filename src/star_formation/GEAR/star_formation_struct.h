/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_GEAR_STAR_FORMATION_STRUCT_H
#define SWIFT_GEAR_STAR_FORMATION_STRUCT_H

/**
 * @brief Functional form of the star formation law
 */
enum star_formation_mode {
  gear_star_formation_default, /*<! Default GEAR star formation mode */
  gear_star_formation_agora    /*<! Agora star formation mode */
};

/* Do we need unique IDs (only useful when spawning
   new particles, conversion gas->stars does not need unique IDs) */
#define star_formation_need_unique_id 1

/**
 * @brief Star-formation-related properties stored in the extended particle
 * data.
 */
struct star_formation_xpart_data {
  /*! Particle velocity divergence. */
  float div_v;
};

/**
 * @brief Star-formation-related properties stored in the star particle
 * data.
 */
struct star_formation_spart_data {

  /*! The birth density */
  float birth_density;

  /*! The birth temperature */
  float birth_temperature;

  /*! The birth mass */
  float birth_mass;

  /*! The progenitor ID */
  long long progenitor_id;
};

/**
 * @brief Global star formation properties
 */
struct star_formation {

  /*! Star formation mode : default or agora */
  int star_formation_mode;

  /*! Number of particle required to resolved the
   * Jeans criterion (at power 2/3). */
  float n_jeans_2_3;

  /*! Maximal gas temperature for forming a star. */
  float maximal_temperature;

  /*! Minimal gas density for forming a star. */
  float density_threshold;

  /*! Star formation efficiency. */
  float star_formation_efficiency;

  /*! Number of possible stars per particle. */
  int n_stars_per_part;

  /*! Mass of a star. */
  float mass_stars;

  /*! Minimal fraction of mass_stars for the last star formed by a part. */
  float min_mass_frac_plus_one;
};

#endif /* SWIFT_GEAR_STAR_FORMATION_STRUCT_H */
