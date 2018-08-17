/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_TRAPHIC_RADIATION_H
#define SWIFT_TRAPHIC_RADIATION_H

/**
 * @brief Initialise the rpart based on the original gas particle.
 *
 * @param p Original gas particle.
 * @param rp #rpart to initialise.
 */
__attribute__((always_inline)) INLINE static void radiation_first_init_part(
    struct part* p, struct rpart* rp) {

  rp->x[0] = p->x[0];
  rp->x[1] = p->x[1];
  rp->x[2] = p->x[2];
  rp->h = p->h;

  rp->density = 1.;
  rp->hydrogen_neutral_fraction = 1.;
  rp->ionising_luminosity = 0.;
}

#endif /* SWIFT_TRAPHIC_RADIATION_H */
