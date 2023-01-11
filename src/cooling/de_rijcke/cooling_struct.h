//
// Created by yuyttenh on 10/01/23.
//

#ifndef SWIFTSIM_COOLING_STRUCT_DE_RIJCKE_H
#define SWIFTSIM_COOLING_STRUCT_DE_RIJCKE_H

/**
 * @brief Properties of the cooling stored in the #part data.
 */
struct cooling_part_data {};

/**
 * @brief Properties of the cooling stored in the extended particle data.
 */
struct cooling_xpart_data {

  /*! Cumulative energy radiated by the particle */
  float radiated_energy;
};

#endif  // SWIFTSIM_COOLING_STRUCT_DE_RIJCKE_H
