/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_STRUCT_GRACKLE_H
#define SWIFT_COOLING_STRUCT_GRACKLE_H

/**
 * @brief Properties of the cooling stored in the #part data.
 */
struct cooling_part_data {};

/**
 * @brief Properties of the cooling stored in the extra particle data
 */
struct cooling_xpart_data {

  /*! Energy radiated away by this particle since the start of the run */
  float radiated_energy;

  /*! Last time the cooling was switch off */
  double time_last_event;

/* here all fractions are mass fraction */
#if COOLING_GRACKLE_MODE >= 1
  /* primordial chemistry >= 1 */
  float HI_frac;
  float HII_frac;
  float HeI_frac;
  float HeII_frac;
  float HeIII_frac;
  float e_frac;

#if COOLING_GRACKLE_MODE >= 2
  /* primordial chemistry >= 2 */
  float HM_frac;
  float H2I_frac;
  float H2II_frac;

#if COOLING_GRACKLE_MODE >= 3
  /* primordial chemistry >= 3 */
  float DI_frac;
  float DII_frac;
  float HDI_frac;
#endif  // MODE >= 3

#endif  // MODE >= 2

#endif  // MODE >= 1

  /* metal cooling = 1 */
  float metal_frac;
};

#endif /* SWIFT_COOLING_STRUCT_GRACKLE_H */
