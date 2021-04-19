/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_COOLING_PROPERTIES_LAMBDA_T_TABLE_H
#define SWIFT_COOLING_PROPERTIES_LAMBDA_T_TABLE_H

/**
 * @file src/cooling/Lambda_T_table/cooling_properties.h
 * @brief Lambda(T) data values.
 */

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! @brief Log T values. */
  double *logT;

  /*! @brief Log Lambda values. */
  double *logLambda;

  /*! @brief Log |dLambda/dT| values. */
  double *logdLambdadT;

  /*! @brief Number of values in the table. */
  int nT;
};

#endif /* SWIFT_COOLING_PROPERTIES_LAMBDA_T_TABLE_H */
