/*******************************************************************************
 * This file is part of GadgetSMP.
 * Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

#ifdef HAVE_HDF5

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#include "task.h"
#include "lock.h"
#include "part.h"
#include "space.h"


/**
 * @brief Reads an HDF5 initial condition file 
 *
 * Opens the HDF5 file fileName and reads the particles contained
 * in the parts array. N is the returned number of particles found
 * in the file.
 *
 */
void read_ic ( char* fileName, struct part *parts,  int* N)
{


}

#endif  /* HAVE_HDF5 */

