
# This file is part of SWIFT.
# Coypright (c) 2012 pedro.gonnet@durham.ac.uk.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# init the project
AC_INIT([SWIFT],[0.1.0])
AC_CONFIG_SRCDIR([src/space.c])
AM_INIT_AUTOMAKE
AC_USE_SYSTEM_EXTENSIONS

LT_INIT

# my own macro collection
AC_CONFIG_MACRO_DIR([m4])

# generate header file
AM_CONFIG_HEADER(config.h)

# find and test the compiler
AC_PROG_CC
AM_PROG_CC_C_O

# Check for MPI. Need to do this before characterising the compiler (C99 mode).
AC_ARG_ENABLE([mpi],
    [AS_HELP_STRING([--enable-mpi],[compile with functionality for distributed-memory parallelim using MPI @<:@default=yes@:>@])],
    [if test "x$enable_mpi" != "xno"
     then
        ACX_MPI([ CC="$MPICC" AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.]) ])
     fi],
    [ACX_MPI([ CC="$MPICC" AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.]) ])])
AM_CONDITIONAL([HAVEMPI],[test -n "$MPICC"])

AC_PROG_CC_C99
AC_C_INLINE

AX_CC_MAXOPT
AX_FUNC_POSIX_MEMALIGN
AX_GCC_ARCHFLAG([no])

#  Check SSE & AVX support.
AX_EXT
if test "$SIMD_FLAGS" != ""; then
   CFLAGS="$CFLAGS $SIMD_FLAGS"
fi

# autoconf stuff
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_HEADER_STDC

# check for the libraries we will need
AC_CHECK_LIB(m,sqrt,,AC_MSG_ERROR(something is wrong with the math library!))

# Check for pthreads
ACX_PTHREAD([LIBS="$PTHREAD_LIBS $LIBS" CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
    CC="$PTHREAD_CC" LDFLAGS="$PTHREAD_LIBS $LIBS"],
    AC_MSG_ERROR([Could not find a working version of
    the pthread library. Make sure you have the library and header files installed
    or use CPPFLAGS and LDFLAGS if the library is installed in a
    non-standard location.]))


# Check for OpenMP
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
if test -z "${OPENMP_CFLAGS}"; then
   echo $OPENMP_CFLAGS
   AC_MSG_ERROR(Compiler does not support OpenMP, 1)
else
   CFLAGS="$CFLAGS $OPENMP_CFLAGS"
fi

# Check for metis
AC_ARG_WITH([metis],
    [AS_HELP_STRING([--with-metis=PATH],[prefix where the metis library is installed @<:@default=yes@:>@])],
    [if test "x$with_metis" != "xno"
     then
        if test "x$with_metis" != "xyes" -a "x$with_metis" != "x"
        then
            METIS_LDFLAGS="-L$with_metis -lmetis"
        else
            METIS_LDFLAGS="-lmetis"
        fi
        AC_CHECK_LIB([metis],[METIS_PartGraphKway],,AC_MSG_ERROR(something is wrong with the metis library!),$METIS_LDFLAGS)
        AC_DEFINE([HAVE_METIS],[true],[The metis library appears to be present.])
	AC_DEFINE(WITH_METIS, 1, [METIS library installed])
     fi])
AC_SUBST(METIS_LDFLAGS)
AM_CONDITIONAL([HAVEMETIS],[test -n "$METIS_LDFLAGS"])


# check for zlib
AC_CHECK_LIB(z,gzopen,[
    AC_DEFINE([HAVE_LIBZ],[1],[Set to 1 if zlib is installed.])
    LDFLAGS="$LDFLAGS -lz"
    ],[])


# Check for HDF5
AX_LIB_HDF5()

# check for setaffinity
AC_CHECK_FUNC( pthread_setaffinity_np , AC_DEFINE([HAVE_SETAFFINITY],[true],
    [Defined if pthread_setaffinity_np exists.]) )

# Check for timing functions needed by cycle.h
AC_HEADER_TIME
AC_CHECK_HEADERS([sys/time.h c_asm.h intrinsics.h mach/mach_time.h])
AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t is defined in <sys/time.h>])],,[#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif])
AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])
AC_MSG_CHECKING([for _rtc intrinsic])
rtc_ok=yes
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif]], [[_rtc()]])],[AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])],[rtc_ok=no])
AC_MSG_RESULT($rtc_ok)

# make the documentation
DX_INIT_DOXYGEN(libswift,doc/Doxyfile,doc/)

# .in files.
AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile doc/Makefile doc/Doxyfile])

# generate output.
AC_OUTPUT
