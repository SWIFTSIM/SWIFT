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

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* This object's header. */
#include "error.h"
#include "version.h"

/* Local headers. */
#include "version_string.h"

/**
 * @brief Return the hostname
 *
 * Will return the name of the host.
 *
 * @result the hostname.
 */
const char *hostname(void) {
  static char buf[256];
  static int initialised = 0;
  if (!initialised) {
    buf[255] = '\0';
    if (gethostname(buf, 255)) sprintf(buf, "%s", "Unknown host");
    initialised = 1;
  }
  return buf;
}

/**
 * @brief Return the source code git revision
 *
 * The SHA of the code checked out when the library was last built.
 * Will include -dirty if they are local modifications.
 *
 * @result the git version.
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
 * The name of the current branch when the code was last built.
 *
 * @result git branch
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
 * @brief Return the date of the commit in the git repository
 *
 * The date of the commit of the code we are running.
 *
 * @result git branch
 */
const char *git_date(void) {
  static char buf[256];
  static int initialised = 0;
  static const char *date = GIT_DATE;
  if (!initialised) {
    if (strlen(date) == 0)
      sprintf(buf, "%s", "unknown");
    else
      sprintf(buf, "%s", date);
    initialised = 1;
  }
  return buf;
}

/**
 * @brief Return the options passed to the 'configure' script
 *
 * @result List of configuration options within simple quotes (').
 */
const char *configuration_options(void) {
  static char buf[1024];
  static int initialised = 0;
  static const char *config = SWIFT_CONFIG_FLAGS;
  if (!initialised) {
    snprintf(buf, 1024, "'%s'", config);
    initialised = 1;
  }
  return buf;
}

/**
 * @brief Return the CFLAGS the code was compiled with
 *
 * @result List of CFLAGS within simple quotes (').
 */
const char *compilation_cflags(void) {
  static char buf[1024];
  static int initialised = 0;
  static const char *cflags = SWIFT_CFLAGS;
  if (!initialised) {
    snprintf(buf, 1024, "'%s'", cflags);
    initialised = 1;
  }
  return buf;
}

/**
 * @brief The version of SWIFT
 *
 * @result the package version
 */
const char *package_version(void) {
  static const char *version = PACKAGE_VERSION;
  return version;
}

/**
 * @brief A description of the package version and code status.
 *
 * @result description of the package version
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

/**
 * @brief return the name of the compiler used to build SWIFT.
 *
 * @result description of the compiler.
 */
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

/**
 * @brief return compiler version used to build SWIFT.
 *
 * @result description of the compiler.
 */
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

/**
 * @brief return the MPI version, runtime if possible otherwise that used when
 *        built.
 *
 * @result description of the MPI version.
 */
const char *mpi_version(void) {
  static char version[80] = {0};

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

  /* Also arbitrarily truncate to keep down to one line, Open MPI,
   * check for last comma and keep to ~60 chars max. */
  strcpy(lib_version + 60, "...");
  ptr = strrchr(lib_version, ',');
  if (ptr != NULL) *ptr = '\0';

#else
  /* Use autoconf guessed value. */
  static char lib_version[60] = {0};
  snprintf(lib_version, 60, "%s", SWIFT_MPI_LIBRARY);
#endif

  /* Numeric version. */
  MPI_Get_version(&std_version, &std_subversion);
  snprintf(version, 80, "%s (MPI std v%i.%i)", lib_version, std_version,
           std_subversion);
#else
  sprintf(version, "Code was not compiled with MPI support");
#endif
  return version;
}

/**
 * @brief return the HDF5 version in use at runtime.
 *
 * @result description of the current HDF5 version.
 */
const char *hdf5_version(void) {

  static char version[256] = {0};
#ifdef HAVE_HDF5
  unsigned int majnum, minnum, relnum;
  H5get_libversion(&majnum, &minnum, &relnum);
  sprintf(version, "%u.%u.%u", majnum, minnum, relnum);
#else
  sprintf(version, "Unknown version");
#endif
  return version;
}

/**
 * @brief return the METIS version used when SWIFT was built.
 *
 * @result description of the METIS version.
 */
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
 * @brief return the FFTW version used when SWIFT was built.
 *
 * @result description of the FFTW version.
 */
const char *fftw3_version(void) {

  static char version[256] = {0};
#if defined(HAVE_FFTW)
  sprintf(version, "%s", "3.x (details not available)");
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

  printf(" Welcome to the cosmological hydrodynamical code\n");
  printf("    ______       _________________\n");
  printf("   / ___/ |     / /  _/ ___/_  __/\n");
  printf("   \\__ \\| | /| / // // /_   / /   \n");
  printf("  ___/ /| |/ |/ // // __/  / /    \n");
  printf(" /____/ |__/|__/___/_/    /_/     \n");
  printf(" SPH With Inter-dependent Fine-grained Tasking\n\n");

  printf(" Version : %s\n", package_version());
  printf(" Revision: %s, Branch: %s, Date: %s\n", git_revision(), git_branch(),
         git_date());
  printf(" Webpage : %s\n\n", PACKAGE_URL);
  printf(" Config. options: %s\n\n", configuration_options());
  printf(" Compiler: %s, Version: %s\n", compiler_name(), compiler_version());
  printf(" CFLAGS  : %s\n", compilation_cflags());
  printf("\n");
#ifdef HAVE_HDF5
  printf(" HDF5 library version: %s\n", hdf5_version());
#endif
#ifdef HAVE_FFTW
  printf(" FFTW library version: %s\n", fftw3_version());
#endif
#ifdef WITH_MPI
  printf(" MPI library: %s\n", mpi_version());
#ifdef HAVE_METIS
  printf(" METIS library version: %s\n", metis_version());
#endif
#endif
  printf("\n");
}
