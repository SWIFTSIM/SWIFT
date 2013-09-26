/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

/* Config parameters. */
#include "../config.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>

#include "const.h"
#include "cycle.h"
#include "part.h"
#include "error.h"
#include "units.h"


/**
 * @brief Initialises the UnitSystem structure with the constants given in const.h
 * @param us The UnitSystem to initialize
 */

void initUnitSystem(struct UnitSystem* us)
{
  us->UnitMass_in_cgs = const_unit_mass_in_cgs;
  us->UnitLength_in_cgs = const_unit_length_in_cgs;
  us->UnitTime_in_cgs = 1.d / (const_unit_velocity_in_cgs / const_unit_length_in_cgs );
  us->UnitCurrent_in_cgs = 1.;
  us->UnitTemperature_in_cgs = 1.;
}


void getBaseUnitExponantsArray(float baseUnitsExp[5], enum UnitConversionFactor unit)
{
  switch( unit )
    {
    case UNIT_CONV_MASS: 
      baseUnitsExp[UNIT_MASS] = 1; break;

    case UNIT_CONV_LENGTH: 
      baseUnitsExp[UNIT_LENGTH] = 1; break;

    case UNIT_CONV_TIME: 
      baseUnitsExp[UNIT_TIME] = 1; break;

    case UNIT_CONV_FREQUENCY: 
       baseUnitsExp[UNIT_TIME] = -1;  break;

    case UNIT_CONV_DENSITY: 
      baseUnitsExp[UNIT_MASS] = 1; baseUnitsExp[UNIT_LENGTH] = -3;  break;

    case UNIT_CONV_SPEED: 
      baseUnitsExp[UNIT_LENGTH] = 1; baseUnitsExp[UNIT_TIME] = -1;  break;

    case UNIT_CONV_ACCELERATION: 
      baseUnitsExp[UNIT_LENGTH] = 1; baseUnitsExp[UNIT_TIME] = -2;  break;

    case UNIT_CONV_FORCE: 
       baseUnitsExp[UNIT_MASS] = 1; baseUnitsExp[UNIT_LENGTH] = 1; baseUnitsExp[UNIT_TIME] = -2;  break;

    case UNIT_CONV_ENERGY: 
       baseUnitsExp[UNIT_MASS] = 1; baseUnitsExp[UNIT_LENGTH] = 2; baseUnitsExp[UNIT_TIME] = -2;  break;

    case UNIT_CONV_ENERGY_PER_UNIT_MASS: 
      baseUnitsExp[UNIT_LENGTH] = 2; baseUnitsExp[UNIT_TIME] = -2;  break;

    case UNIT_CONV_ENTROPY: 
      baseUnitsExp[UNIT_MASS] = 1. - const_hydro_gamma; baseUnitsExp[UNIT_LENGTH] = 3.*const_hydro_gamma - 1.; baseUnitsExp[UNIT_TIME] = -2;  break;

    case UNIT_CONV_POWER: 
      baseUnitsExp[UNIT_MASS] = 1; baseUnitsExp[UNIT_LENGTH] = 2; baseUnitsExp[UNIT_TIME] = -3;  break;

    case UNIT_CONV_PRESSURE: 
      baseUnitsExp[UNIT_MASS] = 1; baseUnitsExp[UNIT_LENGTH] = -1; baseUnitsExp[UNIT_TIME] = -2;  break;

    }
}


/**
 * @brief Returns the conversion factor for a given unit in the chosen unit system
 * @param us The system of units in use
 * @param unit The unit to convert
 */
double conversionFactor(struct UnitSystem* us, enum UnitConversionFactor unit)
{
  float baseUnitsExp[5] = { 0, 0, 0, 0, 0 };

  getBaseUnitExponantsArray(baseUnitsExp, unit);
  
  return generalConversionFactor(us, baseUnitsExp);
}

/**
 * @brief Returns the h factor exponentiation for a given unit
 * @param us The system of units in use
 * @param unit The unit to convert
 */
double hFactor(struct UnitSystem* us, enum UnitConversionFactor unit)
{
  float baseUnitsExp[5] = { 0, 0, 0, 0, 0 };

  getBaseUnitExponantsArray(baseUnitsExp, unit);
  
  return generalhFactor(us, baseUnitsExp);

}


/**
 * @brief Returns the scaling factor exponentiation for a given unit
 * @param us The system of units in use
 * @param unit The unit to convert
 */
double aFactor(struct UnitSystem* us, enum UnitConversionFactor unit)
{
  float baseUnitsExp[5] = { 0, 0, 0, 0, 0 };

  getBaseUnitExponantsArray(baseUnitsExp, unit);
  
  return generalaFactor(us, baseUnitsExp);

}



/**
 * @brief Returns the conversion factor for a given unit (expressed in terms of the 5 fundamental units) in the chosen unit system
 * @param us The unit system used
 * @param baseUnitsExponants The exponant (integer) of each base units required to form the desired quantity. See conversionFactor() for a working example
 */
double generalConversionFactor(struct UnitSystem* us, float baseUnitsExponants[5])
{
  double factor = 1.;
  int i;

  for(i = 0 ; i < 5 ; ++i )
    if(baseUnitsExponants[i] != 0)
      factor *= pow( *( &( us->UnitMass_in_cgs ) + i ) , baseUnitsExponants[i] );

  return factor;	
}


/**
 * @brief Returns the h factor exponentiation for a given unit (expressed in terms of the 5 fundamental units)
 * @param us The unit system used
 * @param baseUnitsExponants The exponant (integer) of each base units required to form the desired quantity. See conversionFactor() for a working example
 */
double generalhFactor(struct UnitSystem* us, float baseUnitsExponants[5])
{
  float factor_exp = 0;
  
  factor_exp += -baseUnitsExponants[UNIT_MASS];
  factor_exp += -baseUnitsExponants[UNIT_LENGTH];
  factor_exp += -baseUnitsExponants[UNIT_TIME];
  

  return factor_exp;	
}

/**
 * @brief Returns the scaling factor exponentiation for a given unit (expressed in terms of the 5 fundamental units)
 * @param us The unit system used
 * @param baseUnitsExponants The exponant (integer) of each base units required to form the desired quantity. See conversionFactor() for a working example
 */
double generalaFactor(struct UnitSystem* us, float baseUnitsExponants[5])
{
  float factor_exp = 0;
  
  factor_exp += baseUnitsExponants[UNIT_LENGTH];
  
  return factor_exp;	
}
