/*******************************************************************************
 * This file is part of SWIFT.
 *
 * The functions in this file are based on code from the HEALPix
 * 3.80 library (see http://healpix.sourceforge.net):
 *
 *  Copyright (C) 1997-2013 Krzysztof M. Gorski, Eric Hivon,
 *                          Benjamin D. Wandelt, Anthony J. Banday,
 *                          Matthias Bartelmann, Hans K. Eriksen,
 *                          Frode K. Hansen, Martin Reinecke
 *
 * Translated and modified for SWIFT by John Helly:
 *
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

#include "lightcone/pixel_index.h"

struct pixel_range {
  pixel_index_t first;
  pixel_index_t last;
};

/*
  Functions we need which are missing from the HEALPix C API
*/

double healpix_max_pixrad(int nside);

void healpix_query_disc_range(int nside, double vec[3], double radius,
                              pixel_index_t *pix_min, pixel_index_t *pix_max,
                              int *nr_ranges, struct pixel_range **range);
