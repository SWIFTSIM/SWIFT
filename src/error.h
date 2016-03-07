/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_ERROR_H
#define SWIFT_ERROR_H

/* Some standard headers. */
#include <stdio.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/**
 * @brief Error macro. Prints the message given in argument and aborts.
 *
 */
#ifdef WITH_MPI
extern int engine_rank;
#define error(s, ...)                                                    \
  {                                                                      \
    fprintf(stderr, "[%03i] %s:%s():%i: " s "\n", engine_rank, __FILE__, \
            __FUNCTION__, __LINE__, ##__VA_ARGS__);                      \
    MPI_Abort(MPI_COMM_WORLD, -1);                                       \
  }
#else
#define error(s, ...)                                                        \
  {                                                                          \
    fprintf(stderr, "%s:%s():%i: " s "\n", __FILE__, __FUNCTION__, __LINE__, \
            ##__VA_ARGS__);                                                  \
    abort();                                                                 \
  }
#endif

#ifdef WITH_MPI
/**
 * @brief MPI error macro. Prints the message given in argument,
 *                         followed by the MPI error string and aborts.
 *
 */
#define mpi_error(res, s, ...)                                           \
  {                                                                      \
    fprintf(stderr, "[%03i] %s:%s():%i: " s "\n", engine_rank, __FILE__, \
            __FUNCTION__, __LINE__, ##__VA_ARGS__);                      \
    int len = 1024;                                                      \
    char buf[len];                                                       \
    MPI_Error_string(res, buf, &len);                                    \
    fprintf(stderr, "%s\n\n", buf);                                      \
    MPI_Abort(MPI_COMM_WORLD, -1);                                       \
  }

#define mpi_error_string(res, s, ...)                                    \
  {                                                                      \
    fprintf(stderr, "[%03i] %s:%s():%i: " s "\n", engine_rank, __FILE__, \
            __FUNCTION__, __LINE__, ##__VA_ARGS__);                      \
    int len = 1024;                                                      \
    char buf[len];                                                       \
    MPI_Error_string(res, buf, &len);                                    \
    fprintf(stderr, "%s\n\n", buf);                                      \
  }
#endif

/**
 * @brief Macro to print a localized message with variable arguments.
 *
 */
#ifdef WITH_MPI
extern int engine_rank;
#define message(s, ...) \
  printf("[%03i] %s: " s "\n", engine_rank, __FUNCTION__, ##__VA_ARGS__)
#else
#define message(s, ...) printf("%s: " s "\n", __FUNCTION__, ##__VA_ARGS__)
#endif

/**
 * @brief Assertion macro compatible with MPI
 *
 */
#ifdef WITH_MPI
extern int engine_rank;
#define assert(expr)                                                       \
  {                                                                        \
    if (!(expr)) {                                                         \
      fprintf(stderr, "[%03i] %s:%s():%i: FAILED ASSERTION: " #expr " \n", \
              engine_rank, __FILE__, __FUNCTION__, __LINE__);              \
      fflush(stderr);                                                      \
      MPI_Abort(MPI_COMM_WORLD, -1);                                       \
    }                                                                      \
  }
#else
#define assert(expr)                                                          \
  {                                                                           \
    if (!(expr)) {                                                            \
      fprintf(stderr, "%s:%s():%i: FAILED ASSERTION: " #expr " \n", __FILE__, \
              __FUNCTION__, __LINE__);                                        \
      fflush(stderr);                                                         \
      abort();                                                                \
    }                                                                         \
  }
#endif

#endif /* SWIFT_ERROR_H */
