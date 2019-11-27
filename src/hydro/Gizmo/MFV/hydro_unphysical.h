/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_HYDRO_UNPHYSICAL_H
#define SWIFT_HYDRO_UNPHYSICAL_H

#if defined(GIZMO_UNPHYSICAL_ERROR) || defined(GIZMO_UNPHYSICAL_RESCUE)

#if defined(GIZMO_UNPHYSICAL_ERROR)

/*! @brief Crash whenever an unphysical value is detected. */
#define gizmo_unphysical_message(name, quantity) \
  error("Unphysical " name " detected (%g)!", quantity);

#elif defined(GIZMO_UNPHYSICAL_WARNING)

/*! @brief Show a warning whenever an unphysical value is detected. */
#define gizmo_unphysical_message(name, quantity) \
  message("Unphysical " name " detected (%g), reset to 0!", quantity);

#else

/*! @brief Don't tell anyone an unphysical value was detected. */
#define gizmo_unphysical_message(name, quantity)

#endif

#define gizmo_check_physical_quantity(name, quantity) \
  if (quantity < 0.f) {                               \
    gizmo_unphysical_message(name, quantity);         \
    quantity = 0.f;                                   \
  }

#define gizmo_check_physical_quantities(                                      \
    mass_name, energy_name, mass, momentum_x, momentum_y, momentum_z, energy) \
  gizmo_check_physical_quantity(mass_name, mass);                             \
  gizmo_check_physical_quantity(energy_name, energy);                         \
  /* now check for vacuum and make sure we have a real vacuum */              \
  if (mass == 0.f || energy == 0.f) {                                         \
    mass = 0.f;                                                               \
    momentum_x = 0.f;                                                         \
    momentum_y = 0.f;                                                         \
    momentum_z = 0.f;                                                         \
    energy = 0.f;                                                             \
  }

#else  // defined(GIZMO_UNPHYSICAL_ERROR) || defined(GIZMO_UNPHYSICAL_RESCUE)

#define gizmo_check_physical_quantities( \
    mass_name, energy_name, mass, momentum_x, momentum_y, momentum_z, energy)

#endif  // defined(GIZMO_UNPHYSICAL_ERROR) || defined(GIZMO_UNPHYSICAL_RESCUE)

#endif /* SWIFT_HYDRO_UNPHYSICAL_H */
