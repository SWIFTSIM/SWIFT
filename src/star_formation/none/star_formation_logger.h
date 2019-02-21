/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_NONE_STAR_FORMATION_LOGGER_H
#define SWIFT_NONE_STAR_FORMATION_LOGGER_H

/* Starformation history struct */
struct star_formation_history {};

INLINE static void starformation_update_SFH(struct spart* sp, struct star_formation_history* sf, const struct cosmology* cosmo, 
    const int with_cosmology){}

INLINE static void starformation_init_SFH(struct star_formation_history* sf, const struct cosmology* cosmo, 
    const int with_cosmology){}

INLINE static void starformation_add_progeny_SFH(struct star_formation_history* sf, 
    const struct star_formation_history* sfprogeny, const struct cosmology* cosmo, 
    const int with_cosmology){}

INLINE static void star_formation_get_total_cell(const struct cell *c, struct star_formation_history *sf, 
    const struct cosmology* cosmo, const int with_cosmology){}


#endif /* SWIFT_NONE_STAR_FORMATION_LOGGER_H */
