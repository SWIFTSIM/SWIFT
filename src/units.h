/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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


/**
 * @brief The unit system used internally.
 *
 * This structure contains the conversion factors to the 7 cgs base units to the internal units.
 * It is used everytime a conversion is performed or an i/o function is called.
 *
 **/
struct UnitSystem
{
  double UnitMass_in_cgs;           /*< Conversion factor from grams to internal mass units */

  double UnitLength_in_cgs;         /*< Conversion factor from centimeters to internal length units. */

  double UnitTime_in_cgs;           /*< Conversion factor from seconds to internal time units. */

  double UnitCurrent_in_cgs;        /*< Conversion factor from Ampere to internal current units. */

  double UnitTemperature_in_cgs;    /*< Conversion factor from Kelvins to internal temperature units. */
};

/**
 * @brief The base units used in the cgs (and internal) system. All units are derived from those.
 */
enum BaseUnits
  {
    UNIT_MASS = 0,
    UNIT_LENGTH = 1,
    UNIT_TIME = 2,
    UNIT_CURRENT = 3,
    UNIT_TEMPERATURE = 4
  };


/**
 * @brief  The different conversion factors supported by default
 */
enum UnitConversionFactor
  {
    UNIT_CONV_NO_UNITS,
    UNIT_CONV_MASS,
    UNIT_CONV_LENGTH,
    UNIT_CONV_TIME,
    UNIT_CONV_DENSITY,
    UNIT_CONV_SPEED,
    UNIT_CONV_ACCELERATION,
    UNIT_CONV_FORCE,
    UNIT_CONV_ENERGY,
    UNIT_CONV_ENERGY_PER_UNIT_MASS,
    UNIT_CONV_ENTROPY,
    UNIT_CONV_ENTROPY_PER_UNIT_MASS,
    UNIT_CONV_POWER,
    UNIT_CONV_PRESSURE,
    UNIT_CONV_FREQUENCY,
    UNIT_CONV_ELECTRIC_CHARGE,
    UNIT_CONV_ELECTRIC_VOLTAGE,
    UNIT_CONV_ELECTRIC_CAPACITANCE,
    UNIT_CONV_ELECTRIC_RESISTANCE,
    UNIT_CONV_ELECTRIC_CONDUCTANCE,
    UNIT_CONV_MAGNETIC_FLUX,
    UNIT_CONV_MAGNETIC_FIELD,
    UNIT_CONV_MAGNETIC_INDUCTANCE,
    UNIT_CONV_TEMPERATURE
  };


/**
 * @brief Initialises the UnitSystem structure with the constants given in const.h
 */
void initUnitSystem(struct UnitSystem*);

/**
 * @brief Returns the base unit conversion factor for a given unit system
 */
double getBaseUnit(struct UnitSystem*, enum BaseUnits);

/**
 * @brief Returns the base unit symbol in the cgs system
 */
const char* getBaseUnitSymbol(enum BaseUnits);

/**
 * @brief Returns the base unit symbol in the cgs system
 */
const char* getBaseUnitCGSSymbol(enum BaseUnits);


/**
 * @brief Returns the conversion factor for a given unit (expressed in terms of the 5 fundamental units) in the chosen unit system
 */
double generalConversionFactor(struct UnitSystem* us, float baseUnitsExponants[5]);


/**
 * @brief Returns the conversion factor for a given unit in the chosen unit system
 */
double conversionFactor(struct UnitSystem* us, enum UnitConversionFactor unit);


/**
 * @brief Returns the h factor for a given unit (expressed in terms of the 5 fundamental units) in the chosen unit system
 */
float generalhFactor(struct UnitSystem* us, float baseUnitsExponants[5]);


/**
 * @brief Returns the h factor for a given unit in the chosen unit system
 */
float hFactor(struct UnitSystem* us, enum UnitConversionFactor unit);


/**
 * @brief Returns the scaling factor for a given unit (expressed in terms of the 5 fundamental units) in the chosen unit system
 */
float generalaFactor(struct UnitSystem* us, float baseUnitsExponants[5]);


/**
 * @brief Returns the scaling factor for a given unit in the chosen unit system
 */
float aFactor(struct UnitSystem* us, enum UnitConversionFactor unit);


/**
 * @brief Returns a string containg the exponants of the base units making up the conversion factors (expressed in terms of the 5 fundamental units)
 */
void generalConversionString(char * buffer, struct UnitSystem* us, float baseUnitsExponants[5]);


/**
 * @brief Returns a string containg the exponants of the base units making up the conversion factors
 */
void conversionString(char * buffer, struct UnitSystem* us, enum UnitConversionFactor unit);
