/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_SOURCETERMS_H
#define SWIFT_SOURCETERMS_H
#include "./const.h"
#ifdef SN_FEEDBACK
#include "sourceterms/sn_feedback/sn_feedback_struct.h"
#endif

/* So far only one model here */
struct sourceterms {
#ifdef SN_FEEDBACK
  struct supernova_struct supernova;
#endif
};

void sourceterms_init(const struct swift_params* parameter_file,
                      struct UnitSystem* us, struct sourceterms* source);
void sourceterms_print(struct sourceterms* source);
#ifdef SN_FEEDBACK
#include "sourceterms/sn_feedback/sn_feedback.h"
#endif
#endif /*  SWIFT_SOURCETERMS_H */
