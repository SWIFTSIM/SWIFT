/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_HYDRO_PACK_H
#define SWIFT_HYDRO_PACK_H

/* Config parameters. */
#include "../config.h"

/* Import the right functions */
#if defined(SPHENIX_SPH)
#include "./hydro/SPHENIX/hydro_pack.h"
#else
#error "No hydro repack functionality for the chosen hydro scheme!"
#endif

#endif /* SWIFT_HYDRO_PACK_H */
