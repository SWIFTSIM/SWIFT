/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_DOIACT_FFT_H
#define SWIFT_RUNNER_DOIACT_FFT_H

struct runner;
struct gravity_tensors;

void runner_do_grav_fft(struct runner* r, int timer);

void multipole_to_mesh_CIC(const struct gravity_tensors* m, double* rho, int N,
                           double fac, const double dim[3]);

void mesh_to_multipole_CIC(struct gravity_tensors* m, const double* pot, int N,
                           double fac, const double dim[3]);

#endif /* SWIFT_RUNNER_DOIACT_FFT_H */
