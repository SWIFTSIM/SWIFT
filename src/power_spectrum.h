/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Marcel van Daalen (daalen@strw.leidenuniv.nl)
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
#ifndef SWIFT_POWER_H
#define SWIFT_POWER_H

/* Config parameters. */
#include <config.h>

/* Local headers */
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

/* Forward declarations */
struct space;
struct gpart;
struct threadpool;
struct swift_params;

/**
 * @brief The different types of components we can calculate the power for.
 */
enum power_type {
  pow_type_matter,
  pow_type_cdm,
  pow_type_gas,
  pow_type_starBH,
  pow_type_pressure,
  pow_type_neutrino,
  pow_type_neutrino_0,
  pow_type_neutrino_1,
  pow_type_count /* Nb. of power spect. types (always last elem. in the enum) */
} __attribute__((packed));

/**
 * @brief Data structure for power spectrum variables and parameters
 */
struct power_spectrum_data {

  /*! Number of grid cells on a side */
  int Ngrid;

  /*! The number of threads used by the FFTW library */
  int nr_threads;

  /*! Number of foldings to use (1 means no folding) */
  int Nfold;

  /*! Factor by which to fold along each dimension */
  int foldfac;

  /*! Number of different power spectra to calculate */
  int spectrumcount;

  /*! The order of the mass assignment window */
  int windoworder;

  /* Shall we correct the position of the k-space bin? */
  int shift_centre_small_k_bins;

  /*! Array of component types to correlate on the "left" side */
  enum power_type* types1;

  /*! Array of component types to correlate on the "right" side */
  enum power_type* types2;

  /*! Pointer to the grid to be reused */
  double* powgrid;

  /*! Pointer to a second grid for cross-power spectra */
  double* powgrid2;

#ifdef HAVE_FFTW
  /*! Pointer to the grid in Fourier space */
  fftw_complex* powgridft;

  /*! Pointer to the second grid in Fourier space */
  fftw_complex* powgridft2;

  /*! The FFT plan to be reused */
  fftw_plan fftplanpow;

  /*! The FFT plan to be reused for the second grid */
  fftw_plan fftplanpow2;
#endif
};

void power_init(struct power_spectrum_data* p, struct swift_params* params,
                int nr_threads);
void calc_all_power_spectra(struct power_spectrum_data* pow_data,
                            const struct space* s, struct threadpool* tp,
                            const int verbose);
void power_clean(struct power_spectrum_data* pow_data);

/* Dump/restore. */
void power_spectrum_struct_dump(const struct power_spectrum_data* p,
                                FILE* stream);
void power_spectrum_struct_restore(struct power_spectrum_data* p, FILE* stream);

#endif /* SWIFT_POWER_H */
