//
// Created by yuyttenh on 25/05/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_UNPHYSICAL_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_UNPHYSICAL_H

#if defined(SHADOWSWIFT_UNPHYSICAL_ERROR) || \
    defined(SHADOWSWIFT_UNPHYSICAL_RESCUE)

#if defined(SHADOWSWIFT_UNPHYSICAL_ERROR)

/*! @brief Crash whenever an unphysical value is detected. */
#define shadowswift_unphysical_message(name, quantity) \
  error("Unphysical " name " detected (%g)!", quantity);

#elif defined(SHADOWSWIFT_UNPHYSICAL_WARNING)

/*! @brief Show a warning whenever an unphysical value is detected. */
#define shadowswift_unphysical_message(name, quantity) \
  warning("Unphysical " name " detected (%g), reset to 0!", quantity);

#else

/*! @brief Don't tell anyone an unphysical value was detected. */
#define shadowswift_unphysical_message(name, quantity)

#endif

#define shadowswift_check_physical_quantity(name, quantity) \
  if (quantity < 0.f) {                                     \
    shadowswift_unphysical_message(name, quantity);         \
    quantity = 0.f;                                         \
  }

#define shadowswift_check_physical_quantities(                          \
    mass_name, energy_name, entropy_name, mass, momentum_x, momentum_y, \
    momentum_z, energy, entropy)                                        \
  shadowswift_check_physical_quantity(mass_name, mass);                 \
  shadowswift_check_physical_quantity(energy_name, energy);             \
  shadowswift_check_physical_quantity(entropy_name, entropy);           \
  /* now check for vacuum and make sure we have a real vacuum */        \
  if (mass == 0.f || energy == 0.f) {                                   \
    mass = 0.f;                                                         \
    momentum_x = 0.f;                                                   \
    momentum_y = 0.f;                                                   \
    momentum_z = 0.f;                                                   \
    energy = 0.f;                                                       \
    entropy = 0.f;                                                      \
  }

#else  // defined(SHADOWSWIFT_UNPHYSICAL_ERROR) ||
       // defined(SHADOWSWIFT_UNPHYSICAL_RESCUE)

#define shadowswift_check_physical_quantities(mass_name, energy_name, mass, \
                                              momentum_x, momentum_y,       \
                                              momentum_z, energy, entropy)

#endif  // defined(SHADOWSWIFT_UNPHYSICAL_ERROR) ||
        // defined(SHADOWSWIFT_UNPHYSICAL_RESCUE)

#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_UNPHYSICAL_H
