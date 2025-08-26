/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DYNAMICAL_FRICTION_STRUCT_BASIC_H
#define SWIFT_DYNAMICAL_FRICTION_STRUCT_BASIC_H


/**
 * @brief DF fields carried by each star particle
 *
 */
struct df_spart_data {

  /* Should DF calculations be done for this particle? */
  /* This will be a mandatory field for all DF models - allows us to skip neighbour finding */
  char apply_df;

  /* Smoothing length for DM neighbours */
  float h_dm;

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density_dm;

  /* Smoothing length for star neighbours */
  float h_stars;

  struct {

    /* Number of neighbours. */
    float wcount;

    /* Number of neighbours spatial derivative. */
    float wcount_dh;

  } density_stars;

  float dm_ngb_mass;

  float dm_rho;

  float dm_v_medium[3];

  float dm_sigma2;

  float dm_df_a[3];

};

#endif /* SWIFT_DYNAMICAL_FRICTION_STRUCT_BASIC_H */
