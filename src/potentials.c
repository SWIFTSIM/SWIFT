#include "units.h"
#include "potentials.h"
#include "physical_constants.h"

void initPotentialProperties(struct UnitSystem *us, struct external_potential* potential)
{
  potential-> point_mass.x    = 50000 * PARSEC_IN_CGS / conversionFactor(us, UNIT_CONV_LENGTH);
  potential-> point_mass.y    = 50000 * PARSEC_IN_CGS / conversionFactor(us, UNIT_CONV_LENGTH);
  potential-> point_mass.z    = 50000 * PARSEC_IN_CGS / conversionFactor(us, UNIT_CONV_LENGTH);
  potential-> point_mass.mass = 1e10 * SOLAR_MASS_IN_CGS / conversionFactor(us, UNIT_CONV_MASS);
}
