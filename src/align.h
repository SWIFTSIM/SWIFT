/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_ALIGN_H
#define SWIFT_ALIGN_H

/**
 * @brief The default struct alignment in SWIFT.
 */
#define SWIFT_STRUCT_ALIGNMENT 32

/**
 * @brief Defines alignment of structures
 */
#define SWIFT_STRUCT_ALIGN __attribute__((aligned(SWIFT_STRUCT_ALIGNMENT)))

/**
 * @brief The default cache alignment in SWIFT.
 */
#define SWIFT_CACHE_ALIGNMENT 64

/**
 * @brief Defines alignment of caches
 */
#define SWIFT_CACHE_ALIGN __attribute__((aligned(SWIFT_CACHE_ALIGNMENT)))

/**
 * @brief Macro to tell the compiler that a given array has the specified
 * alignment.
 *
 * Note that this turns into a no-op but gives information to the compiler.
 *
 * @param array The array.
 * @param alignment The alignment in bytes of the array.
 */
#if defined(__ICC)
#define swift_align_information(array, alignment) \
  __assume_aligned(array, alignment);
#elif defined(__GNUC__)
#define swift_align_information(array, alignment) \
  array = __builtin_assume_aligned(array, alignment);
#else
#define swift_align_information(array, alignment) ;
#endif

/**
 * @brief Macro to create a restrict pointer to an array and tell the compiler
 * that the given array has the specified
 * alignment.
 *
 * Note that this turns into a no-op but gives information to the compiler.
 *
 * @param array The array.
 * @param ptr Pointer to array
 * @param type Type of array
 * @param alignment The alignment in bytes of the array.
 */
#define swift_declare_aligned_ptr(type, array, ptr, alignment) \
  type *restrict array = ptr;                                  \
  swift_align_information(array, alignment);

/**
 * @brief Macro to tell the compiler that a given number is 0 modulo a given
 * size.
 *
 * Note that this turns into a no-op but gives information to the compiler.
 * GCC does not have the equivalent built-in so defaults to nothing.
 *
 * @param var The variable
 * @param size The modulo of interest.
 */
#if defined(__ICC)
#define swift_assume_size(var, size) __assume(var % size == 0);
#else
#define swift_assume_size(var, size) ;
#endif

#endif /* SWIFT_ALIGN_H */
