/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#include "swift.h"

int main() {

  const float h = const_eta_kernel;
  const int numPoints = 30;

  for (int i = 0; i < numPoints; ++i) {

    const float x = i * 3.f / numPoints;
    float W, dW;
    kernel_deval(x / h, &W, &dW);

    printf("h= %f H= %f x=%f W(x,h)=%f\n", h, h * kernel_gamma, x, W);
  }

  return 0;
}
