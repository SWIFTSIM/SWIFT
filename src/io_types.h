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

#include <hdf5.h>

/**
 * @brief Error macro
 *
 */
#define error(s) { fprintf( stderr , "%s:%s():%i: %s\n" , __FILE__ , __FUNCTION__ , __LINE__ , s ); abort(); }

/**
 * @brief The different types of data used in the GADGET IC files.
 *
 * (This is admittedly a poor substitute to C++ templates...)
 */
enum DATA_TYPE{INT, LONG, LONGLONG, UINT, ULONG, ULONGLONG, FLOAT, DOUBLE};

/**
 * @brief The two sorts of data present in the GADGET IC files: compulsory to start a run or optional.
 *
 */
enum DATA_IMPORTANCE{COMPULSORY=1, OPTIONAL=0};

/**
 * @brief Converts a C data type to the HDF5 equivalent. 
 *
 * This function is a trivial wrapper around the HDF5 types but allows
 * to change the exact storage types matching the code types in a transparent way.
 */
static hid_t hdf5Type(enum DATA_TYPE type)
{
  switch(type)
    {
    case INT: return H5T_NATIVE_INT;
    case UINT: return H5T_NATIVE_UINT;
    case LONG: return H5T_NATIVE_LONG;
    case ULONG: return H5T_NATIVE_ULONG;
    case LONGLONG: return H5T_NATIVE_LLONG;
    case ULONGLONG: return H5T_NATIVE_ULLONG;
    case FLOAT: return H5T_NATIVE_FLOAT;
    case DOUBLE: return H5T_NATIVE_DOUBLE;
    default: error("Unknown type");
    }
}

/**
 * @brief Returns the memory size of the data type
 */
static size_t sizeOfType(enum DATA_TYPE type)
{
  switch(type)
    {
    case INT: return sizeof(int);
    case UINT: return sizeof(unsigned int);
    case LONG: return sizeof(long);
    case ULONG: return sizeof(unsigned long);
    case LONGLONG: return sizeof(long long);
    case ULONGLONG: return sizeof(unsigned long long);
    case FLOAT: return sizeof(float);
    case DOUBLE: return sizeof(double);
    default: error("Unknown type");
    }
}


