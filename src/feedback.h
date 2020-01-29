/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_FEEDBACK_H
#define SWIFT_FEEDBACK_H

/* Config parameters. */
#include "../config.h"

/* Select the correct feedback model */
#if defined(FEEDBACK_NONE)
#include "./feedback/none/feedback.h"
#include "./feedback/none/feedback_iact.h"
#elif defined(FEEDBACK_EAGLE)
#include "./feedback/EAGLE/feedback.h"
#include "./feedback/EAGLE/feedback_iact.h"
#elif defined(FEEDBACK_GEAR)
#include "./feedback/GEAR/feedback.h"
#include "./feedback/GEAR/feedback_iact.h"
#else
#error "Invalid choice of feedback model"
#endif

#endif
