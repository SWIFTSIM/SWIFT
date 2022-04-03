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

#include "lightcone/healpix_util.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Integer modulus function, valid for b > 0 only
 *
 * Note that result is not the same as a % b for negative
 * values of a. If b > 0 then result is in range 0 to b-1
 * inclusive.
 *
 * @param a the first integer
 * @param a the second integer
 *
 */
static int mod(int a, int b) {
  int r = a % b;
  return r < 0 ? r + b : r;
}

/**
 * @brief Given a normalized z coordinate, return the ring
 *        number in range 1 to 4*nside-1.
 *
 * @param nside HEALPix resolution parameter
 * @param z the z coordinate to check
 *
 */
static int ring_num(int nside, double z) {

  /* Equatorial regime */
  int iring = round(nside * (2.0 - 1.5 * z));

  /* North cap */
  if (z > 2. / 3.) {
    iring = round(nside * sqrt(3.0 * (1.0 - z)));
    if (iring == 0) iring = 1;
  }

  /* South cap */
  if (z < -2. / 3.) {
    iring = round(nside * sqrt(3.0 * (1.0 + z)));
    if (iring == 0) iring = 1;
    iring = 4 * nside - iring;
  }

  return iring;
}

/**
 * @brief Return information about the specified HEALPix ring
 *
 * @param nside HEALPix resolution parameter
 * @param ring the ring index in range 1 to 4*Nside-1
 * @param npr returns the number of pixels in the ring
 * @param kshift returns the shift of this ring
 * @param npnorth returns total number of pixels in this ring
 *        all rings to the north of this one
 *
 */
static void pixels_per_ring(int nside, int ring, int *npr, int *kshift,
                            long long *npnorth) {

  /* number of pixels in current ring */
  *npr = nside;
  if (ring < *npr) *npr = ring;
  if (4 * nside - ring < *npr) *npr = 4 * nside - ring;
  *npr *= 4;

  /* Shift */
  *kshift = (ring + 1) % 2;              /* 1 for even, 0 for odd */
  if (nside == 1) *kshift = 1 - *kshift; /* except for Nside=1 */
  if (*npr < 4 * nside) *kshift = 1;     /* 1 on polar cap */

  /* Number of pixels in current ring and above */
  if (ring <= nside) {
    /* in North cap */
    *npnorth = ring * (ring + 1ll) * 2ll;
  } else if (ring <= 3 * nside) {
    /* in Equatorial region */
    long long ncap = nside * (nside + 1ll) * 2ll;
    long long ir = ring - nside;
    *npnorth = ncap + 4ll * nside * ir;
  } else {
    /* in South cap */
    long long npix = (12ll * nside) * nside;
    long long ir = 4ll * nside - ring - 1; /* count ring from south */
    *npnorth = npix - ir * (ir + 1ll) * 2ll;
  }
}

/**
 * @brief Compute the z coordinate of a HEALPix ring
 *
 * @param nside HEALPix resolution parameter
 * @param ir the ring index in range 1 to 4*Nside-1
 *
 */
static double ring2z(int nside, int ir) {

  double z;
  double fn = (double)nside;
  if (ir < nside) {
    /* north polar cap */
    double tmp = (double)ir;
    z = 1.0 - (tmp * tmp) / (3.0 * fn * fn);
  } else if (ir < 3 * nside) {
    /* tropical band */
    z = ((double)(2 * nside - ir)) * 2.0 / (3.0 * fn);
  } else {
    /* south polar cap */
    double tmp = (double)(4 * nside - ir);
    z = -1.0 + (tmp * tmp) / (3.0 * fn * fn);
  }
  return z;
}

/**
 * @brief Find pixels with centres within specified radius of
 *        the given vector
 *
 * Based on query_disc() from src/f90/mod/pixel_routines.F90.
 * Assumes RING indexing and does not support inclusive mode
 * (i.e. only returns pixels with centres within radius)
 *
 * If nr_ranges and range are both not NULL, returns a newly
 * allocated array of struct pixel_range with the ranges of
 * pixels which overlap the disc.
 *
 * @param nside HEALPix resolution parameter
 * @param vec vector specifying the disc centre
 * @param radius the radius to search
 * @param pix_min returns minimum pixel index in the disc
 * @param pix_max returns maximum pixel index in the disc
 * @param nr_ranges returns the size of the range array
 * @param range returns a new array of struct pixel_range
 *
 */
void healpix_query_disc_range(int nside, double vec[3], double radius,
                              pixel_index_t *pix_min, pixel_index_t *pix_max,
                              int *nr_ranges, struct pixel_range **range) {

  /* Get normalized disc centre vector */
  double norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  double x0 = vec[0] / norm;
  double y0 = vec[1] / norm;
  double z0 = vec[2] / norm;

  /* coordinate z of highest and lowest points in the disc */
  double rlat0 = asin(z0); /* latitude in RAD of the center */
  double rlat1 = rlat0 + radius;
  double rlat2 = rlat0 - radius;
  double zmin, zmax;
  if (rlat1 >= 0.5 * M_PI) {
    zmax = 1.0;
  } else {
    zmax = sin(rlat1);
  }
  if (rlat2 <= -0.5 * M_PI) {
    zmin = -1.0;
  } else {
    zmin = sin(rlat2);
  }

  /* Find which rings overlap the disc */
  int irmin = ring_num(nside, zmax);
  if (irmin - 1 > 1) {
    irmin = irmin - 1;
  } else {
    irmin = 1;
  }
  int irmax = ring_num(nside, zmin);
  if (irmax + 1 < 4 * nside - 1) {
    irmax = irmax + 1;
  } else {
    irmax = 4 * nside - 1;
  }

  /* Get phi at disc centre */
  double phi0 = 0.0;
  if ((x0 != 0) || (y0 != 0)) phi0 = atan2(y0, x0);

  /* Allocate output array:
     need to allow for worst case where all rings cross the periodic
     boundary and therefore contribute two disjoint ranges of pixels */
  int nout_max = 2 * (irmax - irmin + 1);
  if (nout_max < 1) nout_max = 1;
  if (nr_ranges && range) {
    *range =
        (struct pixel_range *)malloc(nout_max * sizeof(struct pixel_range));
  }

  /* Will return min and max pixel indexes in the disk */
  long long npix = (12ll * nside) * nside;
  long long pix_min_ll = npix;
  long long pix_max_ll = -1;

  /* Now have min/max ring index (in range 1 to 4*nside-1) */
  /* Loop over rings which overlap the disc */
  if (nr_ranges && range) *nr_ranges = 0;
  for (int iring = irmin; iring <= irmax; iring += 1) {

    /* Find z coordinate of this ring */
    double z = ring2z(nside, iring);

    /* Find range in phi which overlaps the disc in this ring:
       taken from  discphirange_at_z() in pix_tools.F90 */
    double cosang = cos(radius);
    double a = x0 * x0 + y0 * y0;
    double dphi = -1000.0; /* Indicates outside disc */
    double b = cosang - z * z0;
    if (a == 0.0) { /* Poles */
      if (b <= 0.0) dphi = M_PI;
    } else {
      double c = fmax(1.0 - z * z, 1.0e-12);
      double cosdphi = b / sqrt(a * c);
      if (cosdphi < -1.0)
        dphi = M_PI; /* all the pixels at this elevation are in the disc */
      if (fabs(cosdphi) <= 1.0) dphi = acos(cosdphi); /* in [0,Pi] */
    }

    /* Look up number of pixels in this ring */
    int npr, kshift;
    long long npnorth;
    pixels_per_ring(nside, iring, &npr, &kshift, &npnorth);

    /* For each ring, store the range of pixels which overlaps the disc.
       If the disc overlaps the periodic boundary at phi=2pi we need to split
       the range into two pieces. */
    int my_low = -1;
    int my_hi = -1;
    if (dphi > M_PI) {
      /* Full ring */
      my_low = 0;
      my_hi = npr - 1;
    } else if (dphi >= 0.0) {
      /* Partial ring */
      double shift = kshift * 0.5;
      int iphi_low = ceil(npr * (phi0 - dphi) / (2 * M_PI) - shift);
      int iphi_hi = floor(npr * (phi0 + dphi) / (2 * M_PI) - shift);
      if (iphi_hi >= iphi_low) {
        my_low = mod(iphi_low, npr);
        my_hi = mod(iphi_hi, npr);
      }
    }
    if (my_low >= 0) {
      long long first;
      long long last;
      if (my_hi >= my_low) {
        /* Not crossing periodic boundary, so we can return a single range */
        first = (npnorth - npr) + my_low;
        last = (npnorth - npr) + my_hi;
        if (first < pix_min_ll) pix_min_ll = first;
        if (last > pix_max_ll) pix_max_ll = last;
        if (nr_ranges && range) {
          (*range)[*nr_ranges].first = first;
          (*range)[*nr_ranges].last = last;
          *nr_ranges += 1;
        }
      } else {
        /* Range overlaps periodic boundary, so will be split in two */
        /* Start of ring to my_hi */
        first = (npnorth - npr) + 0;
        last = (npnorth - npr) + my_hi;
        if (first < pix_min_ll) pix_min_ll = first;
        if (nr_ranges && range) {
          (*range)[*nr_ranges].first = first;
          (*range)[*nr_ranges].last = last;
          *nr_ranges += 1;
        }
        /* my_low to end of ring */
        first = (npnorth - npr) + my_low;
        last = (npnorth - npr) + (npr - 1);
        if (last > pix_max_ll) pix_max_ll = last;
        if (nr_ranges && range) {
          (*range)[*nr_ranges].first = first;
          (*range)[*nr_ranges].last = last;
          *nr_ranges += 1;
        }
      }
    }
    /* Next ring */
  }

  /* Return min and max pixel indexes */
  *pix_min = (pixel_index_t)pix_min_ll;
  *pix_max = (pixel_index_t)pix_max_ll;
}

/**
 * @brief Make a 3D vector given z and phi coordinates
 *
 * @param v returns the new vector
 * @param z normalized coordinate in the z axis
 * @param phi angular coordinate
 *
 */
static void set_z_phi(double *v, double z, double phi) {

  double sintheta = sqrt((1.0 - z) * (1.0 + z));
  v[0] = sintheta * cos(phi);
  v[1] = sintheta * sin(phi);
  v[2] = z;
}

/**
 * @brief Return the maximum radius of any pixel for a given nside.
 *
 * Based on Healpix_base::max_pixrad() from the C++ library.
 *
 * @param nside HEALPix resolution parameter
 *
 */
double healpix_max_pixrad(int nside) {

  double va[3];
  set_z_phi(va, 2. / 3., M_PI / (4 * nside));

  double t1 = 1. - 1. / nside;
  t1 *= t1;
  double vb[3];
  set_z_phi(vb, 1. - t1 / 3., 0.);

  double dotprod = va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2];
  double crossprod[3];
  crossprod[0] = va[1] * vb[2] - va[2] * vb[1];
  crossprod[1] = va[2] * vb[0] - va[0] * vb[2];
  crossprod[2] = va[0] * vb[1] - va[1] * vb[0];
  double length =
      sqrt(crossprod[0] * crossprod[0] + crossprod[1] * crossprod[1] +
           crossprod[2] * crossprod[2]);
  return atan2(length, dotprod);
}
