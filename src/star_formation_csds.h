/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_STAR_FORMATION_PARTICLE_CSDS_H
#define SWIFT_STAR_FORMATION_PARTICLE_CSDS_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "align.h"
#include "csds.h"
#include "part_type.h"
#include "timeline.h"

/* Import the right function */
#if defined(STAR_FORMATION_NONE)
#include "./star_formation/none/star_formation_csds.h"
#elif defined(STAR_FORMATION_QLA)
#error TODO
#elif defined(STAR_FORMATION_EAGLE)
#error TODO
#elif defined(STAR_FORMATION_GEAR)
#include "./star_formation/GEAR/star_formation_csds.h"
#else
#error "Invalid choice of star formation law"
#endif

#endif /* SWIFT_STAR_FORMATION_PARTICLE_CSDS_H */
