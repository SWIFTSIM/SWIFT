/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#include <config.h>

/* Some standard headers. */
#include <fenv.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Conditional headers. */
#ifdef HAVE_LIBZ
#include <zlib.h>
#endif

/* Local headers. */
#include "equation_of_state.h"
#include "swift.h"

/* Engine policy flags. */
#ifndef ENGINE_POLICY
#define ENGINE_POLICY engine_policy_none
#endif

/**
 * @brief Test planetary equations of state.
 *
 *                      WORK IN PROGRESS
 *
 *  *
 */

#ifdef EOS_PLANETARY
int main(int argc, char *argv[]) { return 0; }
#else
int main(int argc, char *argv[]) { return 0; }
#endif
