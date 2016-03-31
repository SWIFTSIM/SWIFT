#ifndef SWIFT_GRAVITY_CONST_H
#define SWIFT_GRAVITY_CONST_H

/* External Potential Constants */

/* Properties of Point Mass */
#ifdef EXTERNAL_POTENTIAL_POINTMASS
#define External_Potential_X  (50000 * PARSEC_IN_CGS /  const_unit_length_in_cgs)
#define External_Potential_Y  (50000 * PARSEC_IN_CGS /  const_unit_length_in_cgs)
#define External_Potential_Z  (50000 * PARSEC_IN_CGS /  const_unit_length_in_cgs)
#define External_Potential_Mass (1e10 * SOLAR_MASS_IN_CGS / const_unit_mass_in_cgs)
#endif

#endif /* SWIFT_GRAVITY_CONST_H */
