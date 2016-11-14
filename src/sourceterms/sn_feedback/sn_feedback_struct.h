/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
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
 * @file src/sourceterms/sn_feedback_struct.h
 * @brief Routines related to source terms (feedback)
 *
 * enumeration type that sets if supernova explosion is done (is_done) or still
 * needs doing (is_not_done)
 */
#ifndef SWIFT_SN_FEEDBACK_STRUCT_H
#define SWIFT_SN_FEEDBACK_STRUCT_H
enum supernova_status { supernova_is_done, supernova_is_not_done };

/**
 * @file src/sourceterms/sn_feedback_struct.h
 * @brief Routines related to source terms (feedback)
 *
 * The structure that describes the source term (supernova feedback)
 * It specifies the time, energy and location of the desired supernova
 * explosion, and a status (supernova_is_done/supernova_is_not_done)
 * that records the status of the supernova
 */
struct supernova_struct {
  double time;
  double energy;
  double x, y, z;
  enum supernova_status status;
};
#endif /* SWIFT_SN_FEEDBACK_STRUCT_H */
