/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl).
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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local headers. */
#include "clocks.h"
#include "memuse.h"

/* Use exit when not developing, avoids core dumps. */
#ifdef SWIFT_DEVELOP_MODE
#define swift_abort(errcode) abort()
#else
#define swift_abort(errcode) exit(errcode)
#endif

/* If reporting memory usage, try to dump that when exiting in error. */
#ifdef SWIFT_MEMUSE_REPORTS
#define memdump(rank) memuse_log_dump_error(rank);
#else
#define memdump(rank)
#endif

/**
 * @brief Error macro. Prints the message given in argument and aborts.
 *
 */
#ifdef WITH_MPI
extern int engine_rank;
#define error(s, ...)                                                      \
  ({                                                                       \
    fflush(stdout);                                                        \
    fprintf(stderr, "[%04i] %s %s:%s():%i: " s "\n", engine_rank,          \
            clocks_get_timesincestart(), __FILE__, __FUNCTION__, __LINE__, \
            ##__VA_ARGS__);                                                \
    memdump(engine_rank);                                                  \
    fflush(stderr);                                                        \
    MPI_Abort(MPI_COMM_WORLD, -1);                                         \
  })
#else
extern int engine_rank;
#define error(s, ...)                                                      \
  ({                                                                       \
    fflush(stdout);                                                        \
    fprintf(stderr, "%s %s:%s():%i: " s "\n", clocks_get_timesincestart(), \
            __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__);              \
    memdump(engine_rank);                                                  \
    fflush(stderr);                                                        \
    swift_abort(1);                                                        \
  })
#endif

#ifdef WITH_MPI
extern int engine_rank;
/**
 * @brief MPI error macro. Prints the message given in argument,
 *                         followed by the MPI error string and aborts.
 *
 */
#define mpi_error(res, s, ...)                                             \
  ({                                                                       \
    fflush(stdout);                                                        \
    fprintf(stderr, "[%04i] %s %s:%s():%i: " s "\n", engine_rank,          \
            clocks_get_timesincestart(), __FILE__, __FUNCTION__, __LINE__, \
            ##__VA_ARGS__);                                                \
    int len = 1024;                                                        \
    char buf[len];                                                         \
    MPI_Error_string(res, buf, &len);                                      \
    fprintf(stderr, "%s\n\n", buf);                                        \
    memdump(engine_rank);                                                  \
    fflush(stderr);                                                        \
    MPI_Abort(MPI_COMM_WORLD, -1);                                         \
  })

#define mpi_error_string(res, s, ...)                                      \
  ({                                                                       \
    fflush(stdout);                                                        \
    fprintf(stderr, "[%04i] %s %s:%s():%i: " s "\n", engine_rank,          \
            clocks_get_timesincestart(), __FILE__, __FUNCTION__, __LINE__, \
            ##__VA_ARGS__);                                                \
    int len = 1024;                                                        \
    char buf[len];                                                         \
    MPI_Error_string(res, buf, &len);                                      \
    fprintf(stderr, "%s\n\n", buf);                                        \
    fflush(stderr);                                                        \
  })
#endif

/**
 * @brief Macro to print a localized message with variable arguments.
 *
 */
#ifdef WITH_MPI
extern int engine_rank;
#define message(s, ...)                                                       \
  ({                                                                          \
    printf("[%04i] %s %s: " s "\n", engine_rank, clocks_get_timesincestart(), \
           __FUNCTION__, ##__VA_ARGS__);                                      \
  })
#else
#define message(s, ...)                                                 \
  ({                                                                    \
    printf("%s %s: " s "\n", clocks_get_timesincestart(), __FUNCTION__, \
           ##__VA_ARGS__);                                              \
  })
#endif

/**
 * @brief Macro to print a simple message with variable arguments.
 * Useful for system messages when time has not started. Assumes these
 * are important, so immediately flushed (can get lost in MPI shutdowns).
 *
 */
#ifdef WITH_MPI
extern int engine_rank;
#define pretime_message(s, ...)                                             \
  ({                                                                        \
    printf("[%04i] %s: " s "\n", engine_rank, __FUNCTION__, ##__VA_ARGS__); \
    fflush(stdout);                                                         \
  })
#else
#define pretime_message(s, ...)                         \
  ({                                                    \
    printf("%s: " s "\n", __FUNCTION__, ##__VA_ARGS__); \
    fflush(stdout);                                     \
  })
#endif

/**
 * @brief Macro to print a localized warning message with variable arguments.
 *
 * Same as message(), but this version prints to the standard error.
 *
 */
#ifdef WITH_MPI
extern int engine_rank;
#define warning(s, ...)                                                \
  ({                                                                   \
    fprintf(stderr, "[%04i] %s %s: WARNING: " s "\n", engine_rank,     \
            clocks_get_timesincestart(), __FUNCTION__, ##__VA_ARGS__); \
  })
#else
#define warning(s, ...)                                                     \
  ({                                                                        \
    fprintf(stderr, "%s %s: WARNING: " s "\n", clocks_get_timesincestart(), \
            __FUNCTION__, ##__VA_ARGS__);                                   \
  })
#endif

/**
 * @brief Assertion macro compatible with MPI
 *
 */
#ifdef WITH_MPI
extern int engine_rank;
#define assert(expr)                                                          \
  ({                                                                          \
    if (!(expr)) {                                                            \
      fflush(stdout);                                                         \
      fprintf(stderr, "[%04i] %s %s:%s():%i: FAILED ASSERTION: " #expr " \n", \
              engine_rank, clocks_get_timesincestart(), __FILE__,             \
              __FUNCTION__, __LINE__);                                        \
      fflush(stderr);                                                         \
      MPI_Abort(MPI_COMM_WORLD, -1);                                          \
    }                                                                         \
  })
#else
#define assert(expr)                                                          \
  ({                                                                          \
    if (!(expr)) {                                                            \
      fflush(stdout);                                                         \
      fprintf(stderr, "%s %s:%s():%i: FAILED ASSERTION: " #expr " \n",        \
              clocks_get_timesincestart(), __FILE__, __FUNCTION__, __LINE__); \
      fflush(stderr);                                                         \
      swift_abort(1);                                                         \
    }                                                                         \
  })
#endif

#ifdef SWIFT_DEBUG_CHECKS

/* Define which cells you'd like to trace. Make them 0 to turn this off. */
#define PROBLEMCELL1 0
#define PROBLEMCELL2 0

/**
 * @brief Macro to trace cells throughout the code.
 */
#ifdef WITH_MPI
extern int engine_rank;
#define celltrace(c, s, ...)                                          \
  ({                                                                  \
    if (c->cellID == PROBLEMCELL1 || c->cellID == PROBLEMCELL2)       \
      printf("[%04i] %s %s: cell %lld local=%d " s "\n", engine_rank, \
             clocks_get_timesincestart(), __FUNCTION__, c->cellID,    \
             c->nodeID == engine_rank, ##__VA_ARGS__);                \
    fflush(stdout);                                                   \
  })
#else
#define celltrace(c, s, ...)                                          \
  ({                                                                  \
    if (c->cellID == PROBLEMCELL1 || c->cellID == PROBLEMCELL2)       \
      printf("%s %s: cell %lld " s "\n", clocks_get_timesincestart(), \
             __FUNCTION__, c->cellID, ##__VA_ARGS__);                 \
    fflush(stdout);                                                   \
  })
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */

/* Define which particles you'd like to trace. */
#define PROBLEMPART1 -1
#define PROBLEMPART2 -1

/**
 * @brief Macro to trace particles throughout the code.
 */
#define parttrace(p, s, ...)                                          \
  ({                                                                  \
    if (p->id == PROBLEMPART1 || p->id == PROBLEMPART2)               \
      printf("%s %s: part %lld " s "\n", clocks_get_timesincestart(), \
             __FUNCTION__, p->id, ##__VA_ARGS__);                     \
    fflush(stdout);                                                   \
  })

#endif /* SWIFT_ERROR_H */
