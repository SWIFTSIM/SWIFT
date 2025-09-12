/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Katy Proctor (katy.proctor@fysik.su.se)
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
#ifndef SWIFT_NONE_SIDM_H
#define SWIFT_NONE_SIDM_H

/**
 * @brief Initialises the si-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sip The particle to act upon
 * @param sidm_properties The properties of the SIDM model.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param scale_factor The current scale-factor (used if running with
 * cosmology).
 * @param time The current time (used if running without cosmology).
 */
__attribute__((always_inline)) INLINE static void sidm_first_init_sipart(
    struct sipart* sip, const struct sidm_props* sidm_properties) {
  sip->time_bin = 0;
}

#endif /* SWIFT_NONE_SIDM_H */
