/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
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

#ifndef SWIFT_DEFAULT_RENDERER_H
#define SWIFT_DEFAULT_RENDERER_H

#include "cosmology.h"
#include "neutrino_properties.h"
#include "physical_constants.h"
#include "units.h"

struct neutrino_renderer {
  double *redshifts;
  double *log_wavenumbers;
  double *ncdm_over_cdm;
  int N_redshifts;
  int N_wavenumbers;
  int N_functions;
};

void renderer_init(struct swift_params *params, const struct unit_system *us,
                   const struct phys_const *phys_const,
                   const struct cosmology *c, const struct neutrino_props *np,
                   struct neutrino_renderer *rend);
void renderer_clean(struct neutrino_renderer *rend);

#endif /* SWIFT_DEFAULT_RENDERER_H */
