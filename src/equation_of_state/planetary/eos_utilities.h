/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_EOS_UTILITIES_H
#define SWIFT_PLANETARY_EOS_UTILITIES_H

/**
 * @file equation_of_state/planetary/eos_utilities.h
 *
 * Utitilies for the planetary equations of state.
 */

/*
    Skip a line while reading a file.
*/
INLINE static int skip_line(FILE *f) {
  int c;

  // Read each character until reaching the end of the line or file
  do {
    c = fgetc(f);
  } while ((c != '\n') && (c != EOF));

  return c;
}

/*
    Skip n lines while reading a file.
*/
INLINE static int skip_lines(FILE *f, int n) {
  int c;

  for (int i = 0; i < n; i++) c = skip_line(f);

  return c;
}

#endif /* SWIFT_PLANETARY_EOS_UTILITIES_H */
