/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Maarten Elion (elion@lorentz.leidenuniv.nl)
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
#ifndef SWIFT_FORCING_IO_H
#define SWIFT_FORCING_IO_H

/**
 * @file src/forcing_struct.h
 * @brief Branches between the different forcing functions 
 */

/* Config parameters. */
#include <config.h>

/* Import the right external potential definition */
#if defined(FORCING_NONE)
#include "./forcing/none/forcing_io.h"
#elif defined(FORCING_ROBERTS_FLOW)
#include "./forcing/roberts_flow/forcing_io.h"
#elif defined(FORCING_ROBERTS_FLOW_ACCELERATION)
#include "./forcing/roberts_flow_acceleration/forcing_io.h"
#elif defined(FORCING_ABC_FLOW)
#include "./forcing/ABC_flow/forcing_io.h"
#elif defined(FORCING_BALSARAKIM)
#include "./forcing/BalsaraKim/forcing_io.h"
#elif defined(FORCING_BOUNDARY_PARTICLES)
#include "./forcing/boundary_particles/forcing_io.h"
#elif defined(FORCING_IDEALIZED_AGN_JET)
#include "./forcing/idealized_agn_jet/forcing_io.h"
#else
#error "Invalid choice of forcing terms"
#endif

#endif /* SWIFT_FORCING_IO_H */