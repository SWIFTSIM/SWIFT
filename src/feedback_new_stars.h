/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FEEDBACK_NEW_STARS_H
#define SWIFT_FEEDBACK_NEW_STARS_H

/* Config parameters. */
#include <config.h>

/* Select the correct feedback model and switches on/off the
 * feedback from stars born in that very step. */
#if defined(FEEDBACK_NONE)
#define feedback_use_newborn_stars 0
#elif defined(FEEDBACK_EAGLE_THERMAL)
#define feedback_use_newborn_stars 0
#elif defined(FEEDBACK_EAGLE_KINETIC)
#define feedback_use_newborn_stars 0
#elif defined(FEEDBACK_GEAR)
#define feedback_use_newborn_stars 1
#elif defined(FEEDBACK_AGORA)
#define feedback_use_newborn_stars 1
#else
#error "Invalid choice of feedback model"
#endif

#endif
