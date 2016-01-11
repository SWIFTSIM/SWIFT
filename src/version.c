/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 * Copyright (C) 2015 Peter W. Draper (p.w.draper@durham.ac.uk).
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

/* Config parameters. */
#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#ifdef HAVE_METIS
#include <metis.h>
#endif
#endif

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

/* Some standard headers. */
#include <stdio.h>
#include <string.h>

/* This object's header. */
#include "version.h"

/**
 * @brief Return the source code git revision
 *
 * @details The SHA of the code checked out when the library was last built.
 * Will include -dirty if they are local modifications.
 */
const char *git_revision(void) {
  static char buf[256];
  static int initialised = 0;
  static const char *revision = GIT_REVISION;
  if (!initialised) {
    if (strlen(revision) == 0)
      sprintf(buf, "%s", "unknown");
    else
      sprintf(buf, "%s", revision);
    initialised = 1;
  }
  return buf;
}

/**
 * @brief Return the source code git branch
 *
 * @details The name of the current branch when the code was last built.
 */
const char *git_branch(void) {
  static char buf[256];
  static int initialised = 0;
  static const char *branch = GIT_BRANCH;
  if (!initialised) {
    if (strlen(branch) == 0)
      sprintf(buf, "%s", "unknown");
    else
      sprintf(buf, "%s", branch);
    initialised = 1;
  }
  return buf;
}

/**
 * @brief The version of SWIFT
 */
const char *package_version(void) {
  static const char *version = PACKAGE_VERSION;
  return version;
}

/**
 * @brief A description of the package version and code status.
 */
const char *package_description(void) {
  static char buf[256];
  static int initialised = 0;
  if (!initialised) {
    sprintf(buf, "SWIFT version: %s, at revision: %s, branch: %s",
            PACKAGE_VERSION, GIT_REVISION, GIT_BRANCH);
    initialised = 1;
  }
  return buf;
}

const char *compiler_name(void) {
  static char compiler[256] = {0};
#if defined(__INTEL_COMPILER)
  sprintf(compiler, "ICC");
#elif defined(__clang__)
  sprintf(compiler, "LLVM/Clang");
#elif defined(__xlc__)
  sprintf(compiler, "IBM XL");
#elif defined(__GNUC__)
  sprintf(compiler, "GCC");
#endif
  if (strlen(compiler) == 0) sprintf(compiler, "Unknown compiler");
  return compiler;
}

const char *compiler_version(void) {
  static char version[256] = {0};
#if defined(__INTEL_COMPILER)
  const int major = __INTEL_COMPILER / 100;
  const int minor = __INTEL_COMPILER - major * 100;
  sprintf(version, "%i.%i.%i", major, minor, __INTEL_COMPILER_BUILD_DATE);
#elif defined(__clang__)
  sprintf(version, "%i.%i.%i", __clang_major__, __clang_minor__,
          __clang_patchlevel__);
#elif defined(__xlc__)
  const int major = __IBMC__ / 100;
  const int minor = (__IBMC__ - major * 100) / 10;
  const int patch = (__IBMC__ - major * 100 - minor * 10);
  sprintf(version, "%i.%i.%i", major, minor, patch);
#elif defined(__GNUC__)
  sprintf(version, "%i.%i.%i", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif
  if (strlen(version) == 0) sprintf(version, "Unknown version");
  return version;
}

const char *mpi_version(void) {
  static char version[256] = {0};
#ifdef WITH_MPI
  int std_version, std_subversion;

  /* Check that the library implements the version string routine */
#ifdef MPI_MAX_LIBRARY_VERSION_STRING
  static char lib_version[MPI_MAX_LIBRARY_VERSION_STRING] = {0};
  int len;
  MPI_Get_library_version(lib_version, &len);

  /* Find first \n and truncate string to this length, can get many lines from
   * some MPIs (MPICH). */
  char *ptr = strchr(lib_version, '\n');
  if (ptr != NULL) *ptr = '\0';

#else
  /* Use autoconf guessed value. */
  static char lib_version[256] = {0};
  sprintf(lib_version, SWIFT_MPI_LIBRARY);
#endif

  /* Numeric version. */
  MPI_Get_version(&std_version, &std_subversion);
  snprintf(version, 256, "%s (implementing MPI standard v %i.%i)", lib_version,
           std_version, std_subversion);
#else
  sprintf(version, "Code was not compiled with MPI support");
#endif
  return version;
}

const char *hdf5_version(void) {

  static char version[256] = {0};
#ifdef HAVE_HDF5
  unsigned int majnum, minnum, relnum;
  H5get_libversion(&majnum, &minnum, &relnum);
  sprintf(version, "%i.%i.%i", majnum, minnum, relnum);
#else
  sprintf(version, "Unknown version");
#endif
  return version;
}

const char *metis_version(void) {

  static char version[256] = {0};
#if defined(WITH_MPI) && defined(HAVE_METIS)
  sprintf(version, "%i.%i.%i", METIS_VER_MAJOR, METIS_VER_MINOR,
          METIS_VER_SUBMINOR);
#else
  sprintf(version, "Unknown version");
#endif
  return version;
}

/**
 * @brief Prints a greeting message to the standard output containing code
 * version and revision number
 */
void greetings(void) {

  printf(" Welcome to the cosmological code\n");
  printf("    ______       _________________\n");
  printf("   / ___/ |     / /  _/ ___/_  __/\n");
  printf("   \\__ \\| | /| / // // /_   / /   \n");
  printf("  ___/ /| |/ |/ // // __/  / /    \n");
  printf(" /____/ |__/|__/___/_/    /_/     \n");
  printf(" SPH With Inter-dependent Fine-grained Tasking\n\n");

  printf(" Version : %s\n", package_version());
  printf(" Revision: %s, Branch: %s\n", git_revision(), git_branch());
  printf(" Webpage : www.swiftsim.com\n\n");
  printf(" Compiler: %s, Version: %s\n", compiler_name(), compiler_version());
#ifdef HAVE_HDF5
  printf(" HDF5 library version: %s\n", hdf5_version());
#endif
#ifdef WITH_MPI
  printf(" MPI library: %s\n", mpi_version());
#ifdef HAVE_METIS
  printf(" METIS library version: %s\n", metis_version());
#endif
#endif
  printf("\n");
}
