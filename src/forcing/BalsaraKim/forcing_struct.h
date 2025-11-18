#ifndef SWIFT_FORCING_BALSARAKIM_STRUCT_H
#define SWIFT_FORCING_BALSARAKIM_STRUCT_H

/**
 * @brief Properties of the cooling stored in the #part data.
 */
struct forcing_part_data {};

/**
 * @brief Properties of the cooling stored in the extended particle data.
 */
struct forcing_xpart_data {

  /*! Cumulative energy injected by forcing */
  float forcing_injected_energy;
};

#endif /* SWIFT_FORCING_BALSARAKIM_STRUCT_H */