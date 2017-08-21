# This file is part of SWIFT.
# Copyright (C) 2012 pedro.gonnet@durham.ac.uk.
#               2016 p.w.draper@durham.ac.uk.
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

# Init the project.
AC_INIT([SWIFT],[0.6.0],[https://gitlab.cosma.dur.ac.uk/swift/swiftsim])
swift_config_flags="$*"

#  Need to define this, instead of using fifth argument of AC_INIT, until 2.64.
AC_DEFINE([PACKAGE_URL],["www.swiftsim.com"], [Package web pages])

AC_COPYRIGHT
AC_CONFIG_SRCDIR([src/space.c])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE

# Add local macro collection.
AC_CONFIG_MACRO_DIR([m4])

# Stop default CFLAGS from anyone except the environment.
: ${CFLAGS=""}

# Generate header file.
AM_CONFIG_HEADER(config.h)

# Find and test the compiler.
AX_CHECK_ENABLE_DEBUG
AC_PROG_CC
AM_PROG_CC_C_O

# Enable POSIX and platform extension preprocessor macros.
AC_USE_SYSTEM_EXTENSIONS

# Check for compiler version and vendor.
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

# Interprocedural optimization support. Needs special handling for linking and
# archiving as well as compilation with Intels, needs to be done before
# libtool is configured (to use correct LD).
AC_ARG_ENABLE([ipo],
   [AS_HELP_STRING([--enable-ipo],
     [Enable interprocedural optimization @<:@no/yes@:>@]
   )],
   [enable_ipo="$enableval"],
   [enable_ipo="no"]
)

if test "$enable_ipo" = "yes"; then
   if test "$ax_cv_c_compiler_vendor" = "intel"; then
      CFLAGS="$CFLAGS -ip -ipo"
      LDFLAGS="$LDFLAGS -ipo"
      : ${AR="xiar"}
      : ${LD="xild"}
      AC_MSG_RESULT([added Intel interprocedural optimization support])
   elif test "$ax_cv_c_compiler_vendor" = "gnu"; then
      CFLAGS="$CFLAGS -flto"
      LDFLAGS="$LDFLAGS -flto"
      AC_MSG_RESULT([added GCC interprocedural optimization support])
   elif test "$ax_cv_c_compiler_vendor" = "clang"; then
      CFLAGS="$CFLAGS -emit-llvm"
      AC_MSG_RESULT([added LLVM interprocedural optimization support])
   else
      AC_MSG_WARN([Compiler does not support interprocedural optimization])
   fi
fi

# Check for MPI. Need to do this before characterising the compiler (C99 mode),
# as this changes the compiler.
# We should consider using AX_PROG_CC_MPI to replace AC_PROG_CC when compiling
# whole applications. There are issues with mixing compilers when using this
# macro. See
# http://lists.gnu.org/archive/html/autoconf-archive-maintainers/2011-05/msg00004.html.
AC_ARG_ENABLE([mpi],
    [AS_HELP_STRING([--enable-mpi],
      [Compile with functionality for distributed-memory parallelism using MPI @<:@yes/no@:>@]
    )],
    [enable_mpi="$enableval"],
    [enable_mpi="yes"]
)
good_mpi="yes"
if test "$enable_mpi" = "yes"; then
    AX_MPI([CC="$MPICC" AC_DEFINE(HAVE_MPI, 1, [Define if you have the MPI library.]) ])
    MPI_LIBRARY="Unknown MPI"

    # Various MPI implementations require additional libraries when also using
    # threads. Use mpirun (on PATH) as that seems to be only command with
    # version flag, allow MPIRUN to override for systems that insist on
    # a non-standard name (PRACE).
    : ${MPIRUN='mpirun'}
    if test "$MPIRUN" = "mpirun"; then
       AC_PATH_PROG([MPIRUN],[mpirun],[notfound])
    fi
    if test "$MPIRUN" = "notfound"; then
       # This may not be fatal (some systems do not allow mpirun on
       # development nodes)., so push on.
       AC_MSG_WARN([Cannot find mpirun command on PATH, thread support may not be correct])
    else
       # Special options we know about.
       # Intel: -mt_mpi
       # PLATFORM: -lmtmpi
       # OpenMPI: nothing, but library should be built correctly.
       # Set MPI_THREAD_LIBS and add to linker commands as necessary.
       AC_MSG_CHECKING([MPI threads options])
       version=`$MPIRUN -version 2>&1`
       case "$version" in
         *Intel*MPI*)
            MPI_THREAD_LIBS="-mt_mpi"
            MPI_LIBRARY="Intel MPI"
            AC_MSG_RESULT([Intel MPI])
         ;;
         *Platform*)
            MPI_THREAD_LIBS="-lmtmpi"
            MPI_LIBRARY="PLATFORM MPI"
            AC_MSG_RESULT([PLATFORM MPI])
         ;;
         *"Open MPI"*)
            MPI_THREAD_LIBS=""
            MPI_LIBRARY="Open MPI"
            AC_MSG_RESULT([Open MPI])
            #  OpenMPI should be 1.8.6 or later, if not complain.
            #  Version is last word on first line of -version output.
            revision=`mpirun -version 2>&1 | grep "Open MPI" | awk '{print $NF}'`
            AX_COMPARE_VERSION( $revision, [ge], [1.8.6],,[good_mpi="no"] )
            if test "$good_mpi" = "no"; then
                AC_MSG_WARN([
    Open MPI version should be at least 1.8.6 (is $revision)])
                enable_mpi="yes (but with warning)"
            fi
         ;;
         *)
            MPI_THREAD_LIBS=""
            AC_MSG_RESULT([unknown])
         ;;
       esac
       AC_SUBST([MPI_THREAD_LIBS])
    fi
    AC_DEFINE_UNQUOTED([SWIFT_MPI_LIBRARY], ["$MPI_LIBRARY"], [The MPI library name, if known.])
fi
AM_CONDITIONAL([HAVEMPI],[test $enable_mpi = "yes"])

# Indicate that MPIRUN can be modified by an environement variable
AC_ARG_VAR(MPIRUN, Path to the mpirun command if non-standard)

# Add libtool support (now that CC is defined).
LT_INIT

# Need C99 and inline support.
AC_PROG_CC_C99
AC_C_INLINE

# If debugging try to show inlined functions.
if test "x$enable_debug" = "xyes"; then
   #  Show inlined functions.
   if test "$ax_cv_c_compiler_vendor" = "gnu"; then
      # Would like to use -gdwarf and let the compiler pick a good version
      # but that doesn't always work.
      AX_CHECK_COMPILE_FLAG([-gdwarf -fvar-tracking-assignments],
        [inline_EXTRA_FLAGS="-gdwarf -fvar-tracking-assignments"],
        [inline_EXTRA_FLAGS="-gdwarf-2 -fvar-tracking-assignments"])
      CFLAGS="$CFLAGS $inline_EXTRA_FLAGS"
   elif test "$ax_cv_c_compiler_vendor" = "intel"; then
      CFLAGS="$CFLAGS -debug inline-debug-info"
   fi
fi

# Check if task debugging is on.
AC_ARG_ENABLE([task-debugging],
   [AS_HELP_STRING([--enable-task-debugging],
     [Store task timing information and generate task dump files @<:@yes/no@:>@]
   )],
   [enable_task_debugging="$enableval"],
   [enable_task_debugging="no"]
)
if test "$enable_task_debugging" = "yes"; then
   AC_DEFINE([SWIFT_DEBUG_TASKS],1,[Enable task debugging])
fi

# Check if threadpool debugging is on.
AC_ARG_ENABLE([threadpool-debugging],
   [AS_HELP_STRING([--enable-threadpool-debugging],
     [Store threadpool mapper timing information and generate threadpool dump files @<:@yes/no@:>@]
   )],
   [enable_threadpool_debugging="$enableval"],
   [enable_threadpool_debugging="no"]
)
if test "$enable_threadpool_debugging" = "yes"; then
   AC_DEFINE([SWIFT_DEBUG_THREADPOOL],1,[Enable threadpool debugging])
   LDFLAGS="$LDFLAGS -rdynamic"
fi

# Check if the general timers are switched on.
AC_ARG_ENABLE([timers],
   [AS_HELP_STRING([--enable-timers],
     [Activate the basic timers @<:@yes/no@:>@]
   )],
   [enable_timers="$enableval"],
   [enable_timers="no"]
)
if test "$enable_timers" = "yes"; then
   AC_DEFINE([SWIFT_USE_TIMERS],1,[Enable individual timers])
fi

# Check if expensive debugging is on.
AC_ARG_ENABLE([debugging-checks],
   [AS_HELP_STRING([--enable-debugging-checks],
     [Activate expensive consistency checks @<:@yes/no@:>@]
   )],
   [enable_debugging_checks="$enableval"],
   [enable_debugging_checks="no"]
)
if test "$enable_debugging_checks" = "yes"; then
   AC_DEFINE([SWIFT_DEBUG_CHECKS],1,[Enable expensive debugging])
fi

# Check if gravity force checks are on for some particles.
AC_ARG_ENABLE([gravity-force-checks],
   [AS_HELP_STRING([--enable-gravity-force-checks],
     [Activate expensive brute-force gravity checks for a fraction 1/N of all particles @<:@N@:>@]
   )],
   [gravity_force_checks="$enableval"],
   [gravity_force_checks="no"]
)
if test "$gravity_force_checks" == "yes"; then
   AC_MSG_ERROR(Need to specify the fraction of particles to check when using --enable-gravity-force-checks!)
elif test "$gravity_force_checks" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_GRAVITY_FORCE_CHECKS], [$enableval] ,[Enable gravity brute-force checks])
fi

# Check if we want to zero the gravity forces for all particles below some ID.
AC_ARG_ENABLE([no-gravity-below-id],
   [AS_HELP_STRING([--enable-no-gravity-below-id],
     [Zeros the gravitational acceleration of all particles with an ID smaller than @<:@N@:>@]
   )],
   [no_gravity_below_id="$enableval"],
   [no_gravity_below_id="no"]
)
if test "$no_gravity_below_id" == "yes"; then
   AC_MSG_ERROR(Need to specify the ID below which particles get zero forces when using --enable-no-gravity-below-id!)
elif test "$no_gravity_below_id" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_NO_GRAVITY_BELOW_ID], [$enableval] ,[Particles with smaller ID than this will have zero gravity forces])
fi


# Define HAVE_POSIX_MEMALIGN if it works.
AX_FUNC_POSIX_MEMALIGN

# Only optimize if allowed, otherwise assume user will set CFLAGS as
# appropriate.
AC_ARG_ENABLE([optimization],
   [AS_HELP_STRING([--enable-optimization],
     [Enable compile time optimization flags for host @<:@yes/no@:>@]
   )],
   [enable_opt="$enableval"],
   [enable_opt="yes"]
)

#  Disable vectorisation for known compilers. This switches off optimizations
#  that could be enabled above, so in general should be appended. Slightly odd
#  implementation as want to describe as --disable-vec, but macro is enable
#  (there is no enable action).
AC_ARG_ENABLE([vec],
   [AS_HELP_STRING([--disable-vec],
     [Disable vectorization]
   )],
   [enable_vec="$enableval"],
   [enable_vec="yes"]
)

HAVEVECTORIZATION=0

if test "$enable_opt" = "yes" ; then

   # Add code optimisation flags and tuning to host. This is a funny macro
   # that does not like CFLAGS being already set. Work around that as we have
   # at least set it to "", so it is set.
   ac_test_CFLAGS="no"
   old_CFLAGS="$CFLAGS"
   AX_CC_MAXOPT
   ac_test_CFLAGS="yes"
   CFLAGS="$old_CFLAGS $CFLAGS"

   # Check SSE & AVX support (some overlap with AX_CC_MAXOPT).
   # Don't use the SIMD_FLAGS result with Intel compilers. The -x<code>
   # value from AX_CC_MAXOPT should be sufficient.
   AX_EXT
   if test "$SIMD_FLAGS" != ""; then
       if test "$ax_cv_c_compiler_vendor" != "intel"; then
           CFLAGS="$CFLAGS $SIMD_FLAGS"
       fi
   fi


   if test "$enable_vec" = "no"; then
      if test "$ax_cv_c_compiler_vendor" = "intel"; then
      	 CFLAGS="$CFLAGS -no-vec -no-simd"
      	 AC_MSG_RESULT([disabled Intel vectorization])
      elif test "$ax_cv_c_compiler_vendor" = "gnu"; then
      	 CFLAGS="$CFLAGS -fno-tree-vectorize"
      	 AC_MSG_RESULT([disabled GCC vectorization])
      elif test "$ax_cv_c_compiler_vendor" = "clang"; then
         CFLAGS="$CFLAGS -fno-vectorize -fno-slp-vectorize"
         AC_MSG_RESULT([disabled clang vectorization])
      else
         AC_MSG_WARN([Do not know how to disable vectorization for this compiler])
      fi
   else
      AC_DEFINE([WITH_VECTORIZATION],1,[Enable vectorization])
      HAVEVECTORIZATION=1
   fi
fi
AM_CONDITIONAL([HAVEVECTORIZATION],[test -n "$HAVEVECTORIZATION"])


# Add address sanitizer options to flags, if requested. Only useful for GCC
# version 4.8 and later and clang.
AC_ARG_ENABLE([sanitizer],
   [AS_HELP_STRING([--enable-sanitizer],
     [Enable memory error detection using address sanitizer @<:@no/yes@:>@]
   )],
   [enable_san="$enableval"],
   [enable_san="no"]
)

if test "$enable_san" = "yes"; then
   if test "$ax_cv_c_compiler_vendor" = "gnu"; then
      AX_COMPARE_VERSION( $ax_cv_c_compiler_version, [ge], [4.8.0],
                          [enable_san="yes"], [enable_san="no"] )
   elif test "$ax_cv_c_compiler_vendor" = "clang"; then
      AX_COMPARE_VERSION( $ax_cv_c_compiler_version, [ge], [3.2.0],
                          [enable_san="yes"], [enable_san="no"] )
   fi
   if test "$enable_san" = "yes"; then
      CFLAGS="$CFLAGS -fsanitize=address -fno-omit-frame-pointer"
      AC_MSG_RESULT([added address sanitizer support])
   else
      AC_MSG_WARN([Compiler does not support address sanitizer option])
   fi
fi

# Autoconf stuff.
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_HEADER_STDC

# Check for the libraries we will need.
AC_CHECK_LIB(m,sqrt,,AC_MSG_ERROR(something is wrong with the math library!))

# Check for pthreads.
AX_PTHREAD([LIBS="$PTHREAD_LIBS $LIBS" CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
    CC="$PTHREAD_CC" LDFLAGS="$LDFLAGS $PTHREAD_LIBS $LIBS"],
    AC_MSG_ERROR([Could not find a working version of
    the pthread library. Make sure you have the library and header files installed
    or use CPPFLAGS and LDFLAGS if the library is installed in a
    non-standard location.]))


# Check for metis. Note AX_LIB_METIS exists, but cannot be configured
# to be default off (i.e. given no option it tries to locate METIS), so we
# don't use that.
have_metis="no"
AC_ARG_WITH([metis],
    [AS_HELP_STRING([--with-metis=PATH],
       [root directory where metis is installed @<:@yes/no@:>@]
    )],
    [],
    [with_metis="no"]
)
if test "x$with_metis" != "xno"; then
   if test "x$with_metis" != "xyes" -a "x$with_metis" != "x"; then
      METIS_LIBS="-L$with_metis/lib -lmetis"
      METIS_INCS="-I$with_metis/include"
   else
      METIS_LIBS="-lmetis"
      METIS_INCS=""
   fi
   have_metis="yes"
   AC_CHECK_LIB([metis],[METIS_PartGraphKway],
      AC_DEFINE([HAVE_METIS],1,[The metis library appears to be present.]),
      AC_MSG_ERROR(something is wrong with the metis library!),$METIS_LIBS)
fi
AC_SUBST([METIS_LIBS])
AC_SUBST([METIS_INCS])
AM_CONDITIONAL([HAVEMETIS],[test -n "$METIS_LIBS"])

#  Check for tcmalloc a fast malloc that is part of the gperftools.
have_tcmalloc="no"
AC_ARG_WITH([tcmalloc],
   [AS_HELP_STRING([--with-tcmalloc],
      [use tcmalloc library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_tcmalloc="$withval"],
   [with_tcmalloc="no"]
)
if test "x$with_tcmalloc" != "xno"; then
   if test "x$with_tcmalloc" != "xyes" -a "x$with_tcmalloc" != "x"; then
      tclibs="-L$with_tcmalloc -ltcmalloc"
   else
      tclibs="-ltcmalloc"
   fi
   AC_CHECK_LIB([tcmalloc],[tc_cfree],[have_tcmalloc="yes"],[have_tcmalloc="no"],
                $tclibs)

   #  Could just have the minimal version.
   if test "$have_tcmalloc" = "no"; then
      if test "x$with_tcmalloc" != "xyes" -a "x$with_tcmalloc" != "x"; then
         tclibs="-L$with_tcmalloc -ltcmalloc_minimal"
      else
         tclibs="-ltcmalloc_minimal"
      fi
      AC_CHECK_LIB([tcmalloc],[tc_cfree],[have_tcmalloc="yes"],[have_tcmalloc="no"],
                   $tclibs)
   fi

   if test "$have_tcmalloc" = "yes"; then
      TCMALLOC_LIBS="$tclibs"

      # These are recommended for GCC.
      if test "$ax_cv_c_compiler_vendor" = "gnu"; then
         CFLAGS="$CFLAGS -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free"
      fi
   else
      TCMALLOC_LIBS=""
   fi
fi
AC_SUBST([TCMALLOC_LIBS])
AM_CONDITIONAL([HAVETCMALLOC],[test -n "$TCMALLOC_LIBS"])

#  Check for -lprofiler usually part of the gperftools along with tcmalloc.
have_profiler="no"
AC_ARG_WITH([profiler],
   [AS_HELP_STRING([--with-profiler],
      [use cpu profiler library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_profiler="$withval"],
   [with_profiler="yes"]
)
if test "x$with_profiler" != "xno"; then
   if test "x$with_profiler" != "xyes" -a "x$with_profiler" != "x"; then
      proflibs="-L$with_profiler -lprofiler"
   else
      proflibs="-lprofiler"
   fi
   AC_CHECK_LIB([profiler],[ProfilerFlush],[have_profiler="yes"],[have_profiler="no"],
                $proflibs)

   if test "$have_profiler" = "yes"; then
      PROFILER_LIBS="$proflibs"
   else
      PROFILER_LIBS=""
   fi
fi
AC_SUBST([PROFILER_LIBS])
AM_CONDITIONAL([HAVEPROFILER],[test -n "$PROFILER_LIBS"])

#  Check for jemalloc another fast malloc that is good with contention.
have_jemalloc="no"
AC_ARG_WITH([jemalloc],
   [AS_HELP_STRING([--with-jemalloc],
      [use jemalloc library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_jemalloc="$withval"],
   [with_jemalloc="no"]
)
if test "x$with_jemalloc" != "xno"; then
   if test "x$with_jemalloc" != "xyes" -a "x$with_jemalloc" != "x"; then
      jelibs="-L$with_jemalloc -ljemalloc"
   else
      jelibs="-ljemalloc"
   fi
   AC_CHECK_LIB([jemalloc],[malloc_usable_size],[have_jemalloc="yes"],[have_jemalloc="no"],
                $jelibs)

   if test "$have_jemalloc" = "yes"; then
      JEMALLOC_LIBS="$jelibs"
   else
      JEMALLOC_LIBS=""
   fi
fi
AC_SUBST([JEMALLOC_LIBS])
AM_CONDITIONAL([HAVEJEMALLOC],[test -n "$JEMALLOC_LIBS"])

#  Don't allow both tcmalloc and jemalloc.
if test "x$have_tcmalloc" != "xno" -a "x$have_jemalloc" != "xno"; then
   AC_MSG_ERROR([Cannot use tcmalloc at same time as jemalloc])
fi

# Check for HDF5. This is required.
AX_LIB_HDF5
if test "$with_hdf5" != "yes"; then
    AC_MSG_ERROR([Could not find a working HDF5 library])
fi

# We want to know if this HDF5 supports MPI and whether we should use it.
# The default is to use MPI support if it is available, i.e. this is
# a parallel HDF5.
have_parallel_hdf5="no"
if test "$with_hdf5" = "yes"; then
    AC_ARG_ENABLE([parallel-hdf5],
       [AS_HELP_STRING([--enable-parallel-hdf5],
         [Enable parallel HDF5 library MPI functions if available. @<:@yes/no@:>@]
       )],
       [enable_parallel_hdf5="$enableval"],
       [enable_parallel_hdf5="yes"]
    )

    if test "$enable_parallel_hdf5" = "yes"; then
        AC_MSG_CHECKING([for HDF5 parallel support])
# Check if the library is capable, the header should define H5_HAVE_PARALLEL.

        AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
        #include "hdf5.h"
        #ifndef H5_HAVE_PARALLEL
        # error macro not defined
        #endif
        ]])], [parallel="yes"], [parallel="no"])
        if test "$parallel" = "yes"; then
            have_parallel_hdf5="yes"
            AC_DEFINE([HAVE_PARALLEL_HDF5],1,[HDF5 library supports parallel access])
        fi
        AC_MSG_RESULT($parallel)
    fi
fi
AM_CONDITIONAL([HAVEPARALLELHDF5],[test "$have_parallel_hdf5" = "yes"])

# Check for setaffinity.
AC_CHECK_FUNC(pthread_setaffinity_np, AC_DEFINE([HAVE_SETAFFINITY],[1],
    [Defined if pthread_setaffinity_np exists.]) )
AM_CONDITIONAL(HAVESETAFFINITY,
    [test "$ac_cv_func_pthread_setaffinity_np" = "yes"])

have_numa="no"
if test "$ac_cv_func_pthread_setaffinity_np" = "yes"; then
  # Check for libnuma.
  AC_CHECK_HEADER([numa.h])
  if test "$ac_cv_header_numa_h" = "yes"; then
    AC_CHECK_LIB([numa], [numa_available])
    have_numa="yes"
  fi
fi

# Check for FFTW
have_fftw3="no"
AC_CHECK_HEADER([fftw3.h])
if test "$ac_cv_header_fftw3_h" = "yes"; then
   AC_CHECK_LIB([fftw3],[fftw_malloc], [AC_DEFINE([HAVE_FFTW],
   [1],[Defined if FFTW 3.x exists.])] )
   FFTW_LIBS="-lfftw3"
   have_fftw3="yes"
fi
AC_SUBST([FFTW_LIBS])

# Check for Intel and PowerPC intrinsics header optionally used by vector.h.
AC_CHECK_HEADERS([immintrin.h])
AC_CHECK_HEADERS([altivec.h])

# Check for timing functions needed by cycle.h.
AC_HEADER_TIME
AC_CHECK_HEADERS([sys/time.h c_asm.h intrinsics.h mach/mach_time.h])
AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t
is defined in <sys/time.h>])],,
[#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif])
AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])
AC_MSG_CHECKING([for _rtc intrinsic])
rtc_ok=yes
AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif]],
[[_rtc()]])],
[AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])],[rtc_ok=no])
AC_MSG_RESULT($rtc_ok)

# Add warning flags by default, if these can be used. Option =error adds
# -Werror to GCC, clang and Intel.  Note do this last as compiler tests may
# become errors, if that's an issue don't use CFLAGS for these, use an AC_SUBST().
AC_ARG_ENABLE([compiler-warnings],
   [AS_HELP_STRING([--enable-compiler-warnings],
     [Enable compile time warning flags, if compiler is known @<:@error/no/yes)@:>@]
   )],
   [enable_warn="$enableval"],
   [enable_warn="error"]
)
if test "$enable_warn" != "no"; then

    # AX_CFLAGS_WARN_ALL does not give good warning flags for the Intel compiler
    # We will do this by hand instead and only default to the macro for unknown compilers
    case "$ax_cv_c_compiler_vendor" in
          gnu | clang)
             CFLAGS="$CFLAGS -Wall -Wextra -Wno-unused-parameter"
          ;;
	  intel)
             CFLAGS="$CFLAGS -w2 -Wunused-variable"
          ;;
	  *)
	     AX_CFLAGS_WARN_ALL
	  ;;
    esac

    # Add a "choke on warning" flag if it exists
    if test "$enable_warn" = "error"; then
       case "$ax_cv_c_compiler_vendor" in
          intel | gnu | clang)
             CFLAGS="$CFLAGS -Werror"
          ;;
       esac
    fi
fi

# Various package configuration options.

# Hydro scheme.
AC_ARG_WITH([hydro],
   [AS_HELP_STRING([--with-hydro=<scheme>],
      [Hydro dynamics to use @<:@gadget2, minimal, hopkins, default, gizmo, shadowfax default: gadget2@:>@]
   )],
   [with_hydro="$withval"],
   [with_hydro="gadget2"]
)
case "$with_hydro" in
   gadget2)
      AC_DEFINE([GADGET2_SPH], [1], [Gadget-2 SPH])
   ;; 
   minimal)
      AC_DEFINE([MINIMAL_SPH], [1], [Minimal SPH])
   ;; 
   hopkins)
      AC_DEFINE([HOPKINS_PE_SPH], [1], [Pressure-Entropy SPH])
   ;; 
   default)
      AC_DEFINE([DEFAULT_SPH], [1], [Default SPH])
   ;; 
   gizmo)
      AC_DEFINE([GIZMO_SPH], [1], [GIZMO SPH])
   ;; 
   shadowfax)
      AC_DEFINE([SHADOWFAX_SPH], [1], [Shadowfax SPH])
   ;; 

   *)
      AC_MSG_ERROR([Unknown hydrodynamics scheme: $with_hydro])
   ;;
esac

# SPH Kernel function
AC_ARG_WITH([kernel],
   [AS_HELP_STRING([--with-kernel=<kernel>],
      [Kernel function to use @<:@cubic-spline, quartic-spline, quintic-spline, wendland-C2, wendland-C4, wendland-C6 default: cubic-spline@:>@]
   )],
   [with_kernel="$withval"],
   [with_kernel="cubic-spline"]
)
case "$with_kernel" in
   cubic-spline)
      AC_DEFINE([CUBIC_SPLINE_KERNEL], [1], [Cubic spline kernel])
   ;;    
   quartic-spline)
      AC_DEFINE([QUARTIC_SPLINE_KERNEL], [1], [Quartic spline kernel])
   ;;    
   quintic-spline)
      AC_DEFINE([QUINTIC_SPLINE_KERNEL], [1], [Quintic spline kernel])
   ;;    
   wendland-C2)
      AC_DEFINE([WENDLAND_C2_KERNEL], [1], [Wendland-C2 kernel])
   ;;    
   wendland-C4)
      AC_DEFINE([WENDLAND_C4_KERNEL], [1], [Wendland-C4 kernel])
   ;;    
   wendland-C6)
      AC_DEFINE([WENDLAND_C6_KERNEL], [1], [Wendland-C6 kernel])
   ;;    
   *)
      AC_MSG_ERROR([Unknown kernel function: $with_kernel])
   ;;
esac

#  Dimensionality of the hydro scheme.
AC_ARG_WITH([hydro-dimension],
   [AS_HELP_STRING([--with-hydro-dimension=<dim>],
      [dimensionality of problem @<:@3/2/1 default: 3@:>@]
   )],
   [with_dimension="$withval"],
   [with_dimension="3"]
)
case "$with_dimension" in
   1)
      AC_DEFINE([HYDRO_DIMENSION_1D], [1], [1D solver])
   ;; 
   2)
      AC_DEFINE([HYDRO_DIMENSION_2D], [2], [2D solver])
   ;; 
   3)
      AC_DEFINE([HYDRO_DIMENSION_3D], [3], [3D solver])
   ;; 
   *)
      AC_MSG_ERROR([Dimensionality must be 1, 2 or 3])
   ;;
esac

#  Equation of state
AC_ARG_WITH([equation-of-state],
   [AS_HELP_STRING([--with-equation-of-state=<EoS>],
      [equation of state @<:@ideal-gas, isothermal-gas default: ideal-gas@:>@]
   )],
   [with_eos="$withval"],
   [with_eos="ideal-gas"]
)
case "$with_eos" in
   ideal-gas)
      AC_DEFINE([EOS_IDEAL_GAS], [1], [Ideal gas equation of state])
   ;; 
   isothermal-gas)
      AC_DEFINE([EOS_ISOTHERMAL_GAS], [1], [Isothermal gas equation of state])
   ;; 
   *)
      AC_MSG_ERROR([Unknown equation of state: $with_eos])
   ;;
esac

#  Adiabatic index
AC_ARG_WITH([adiabatic-index],
   [AS_HELP_STRING([--with-adiabatic-index=<gamma>],
      [adiabatic index @<:@5/3, 7/5, 4/3, 2 default: 5/3@:>@]
   )],
   [with_gamma="$withval"],
   [with_gamma="5/3"]
)
case "$with_gamma" in
   5/3)
      AC_DEFINE([HYDRO_GAMMA_5_3], [5./3.], [Adiabatic index is 5/3])
   ;; 
   7/5)
      AC_DEFINE([HYDRO_GAMMA_7_5], [7./5.], [Adiabatic index is 7/5])
   ;; 
   4/3)
      AC_DEFINE([HYDRO_GAMMA_4_3], [4./3.], [Adiabatic index is 4/3])
   ;; 
   2)
      AC_DEFINE([HYDRO_GAMMA_2_1], [2.], [Adiabatic index is 2])
   ;; 
   *)
      AC_MSG_ERROR([Unknown adiabatic index: $with_gamma])
   ;;
esac

#  Riemann solver
AC_ARG_WITH([riemann-solver],
   [AS_HELP_STRING([--with-riemann-solver=<solver>],
      [riemann solver (gizmo-sph only) @<:@none, exact, trrs, hllc, default: none@:>@]
   )],
   [with_riemann="$withval"],
   [with_riemann="none"]
)
case "$with_riemann" in
   none)
      AC_DEFINE([RIEMANN_SOLVER_NONE], [1], [No Riemann solver])
   ;; 
   exact)
      AC_DEFINE([RIEMANN_SOLVER_EXACT], [1], [Exact Riemann solver])
   ;; 
   trrs)
      AC_DEFINE([RIEMANN_SOLVER_TRRS], [1], [Two Rarefaction Riemann Solver])
   ;; 
   hllc)
      AC_DEFINE([RIEMANN_SOLVER_HLLC], [1], [Harten-Lax-van Leer-Contact Riemann solver])
   ;; 
   *)
      AC_MSG_ERROR([Unknown Riemann solver: $with_riemann])
   ;;
esac

#  Cooling function
AC_ARG_WITH([cooling],
   [AS_HELP_STRING([--with-cooling=<function>],
      [cooling function @<:@none, const-du, const-lambda, grackle default: none@:>@]
   )],
   [with_cooling="$withval"],
   [with_cooling="none"]
)
case "$with_cooling" in
   none)
      AC_DEFINE([COOLING_NONE], [1], [No cooling function])
   ;; 
   const-du)
      AC_DEFINE([COOLING_CONST_DU], [1], [Const du/dt cooling function])
   ;; 
   const-lambda)
      AC_DEFINE([COOLING_CONST_LAMBDA], [1], [Const Lambda cooling function])
   ;; 
   grackle)
      AC_DEFINE([COOLING_GRACKLE], [1], [Cooling via the grackle library])
   ;; 
   *)
      AC_MSG_ERROR([Unknown cooling function: $with_cooling])
   ;;
esac

#  External potential
AC_ARG_WITH([ext-potential],
   [AS_HELP_STRING([--with-ext-potential=<pot>],
      [external potential @<:@none, point-mass, isothermal, softened-isothermal, disc-patch, sine-wave default: none@:>@]
   )],
   [with_potential="$withval"],
   [with_potential="none"]
)
case "$with_potential" in
   none)
      AC_DEFINE([EXTERNAL_POTENTIAL_NONE], [1], [No external potential])
   ;; 
   point-mass)
      AC_DEFINE([EXTERNAL_POTENTIAL_POINTMASS], [1], [Point-mass external potential])
   ;; 
   isothermal)
      AC_DEFINE([EXTERNAL_POTENTIAL_ISOTHERMAL], [1], [Isothermal external potential])
   ;; 
   disc-patch)
      AC_DEFINE([EXTERNAL_POTENTIAL_DISC_PATCH], [1], [Disc-patch external potential])
   ;; 
   sine-wave)
      AC_DEFINE([EXTERNAL_POTENTIAL_SINE_WAVE], [1], [Sine wave external potential in 1D])
   ;; 
   *)
      AC_MSG_ERROR([Unknown external potential: $with_potential])
   ;;
esac

#  Gravity multipole order
AC_ARG_WITH([multipole-order],
   [AS_HELP_STRING([--with-multipole-order=<order>],
      [order of the multipole and gravitational field expansion @<:@ default: 4@:>@]
   )],
   [with_multipole_order="$withval"],
   [with_multipole_order="4"]
)
AC_DEFINE_UNQUOTED([SELF_GRAVITY_MULTIPOLE_ORDER], [$with_multipole_order], [Multipole order])


# Check for git, needed for revision stamps.
AC_PATH_PROG([GIT_CMD], [git])
AC_SUBST([GIT_CMD])

# Make the documentation. Add conditional to handle disable option.
DX_INIT_DOXYGEN(libswift,doc/Doxyfile,doc/)
AM_CONDITIONAL([HAVE_DOXYGEN], [test "$ac_cv_path_ac_pt_DX_DOXYGEN" != ""])

# Handle .in files.
AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile doc/Makefile doc/Doxyfile tests/Makefile])
AC_CONFIG_FILES([tests/testReading.sh], [chmod +x tests/testReading.sh])
AC_CONFIG_FILES([tests/testActivePair.sh], [chmod +x tests/testActivePair.sh])
AC_CONFIG_FILES([tests/test27cells.sh], [chmod +x tests/test27cells.sh])
AC_CONFIG_FILES([tests/test27cellsPerturbed.sh], [chmod +x tests/test27cellsPerturbed.sh])
AC_CONFIG_FILES([tests/test125cells.sh], [chmod +x tests/test125cells.sh])
AC_CONFIG_FILES([tests/test125cellsPerturbed.sh], [chmod +x tests/test125cellsPerturbed.sh])
AC_CONFIG_FILES([tests/testPeriodicBC.sh], [chmod +x tests/testPeriodicBC.sh])
AC_CONFIG_FILES([tests/testPeriodicBCPerturbed.sh], [chmod +x tests/testPeriodicBCPerturbed.sh])
AC_CONFIG_FILES([tests/testInteractions.sh], [chmod +x tests/testInteractions.sh])
AC_CONFIG_FILES([tests/testParser.sh], [chmod +x tests/testParser.sh])

# Save the compilation options
AC_DEFINE_UNQUOTED([SWIFT_CONFIG_FLAGS],["$swift_config_flags"],[Flags passed to configure])

# Make sure the latest git revision string gets included
touch src/version.c

# Generate output.
AC_OUTPUT

# Report general configuration.
AC_MSG_RESULT([ 
 ------- Summary --------

   $PACKAGE_NAME v.$PACKAGE_VERSION

   Compiler        : $CC
    - vendor       : $ax_cv_c_compiler_vendor
    - version      : $ax_cv_c_compiler_version
    - flags        : $CFLAGS
   MPI enabled     : $enable_mpi
   HDF5 enabled    : $with_hdf5
    - parallel     : $have_parallel_hdf5
   Metis enabled   : $have_metis
   FFTW3 enabled   : $have_fftw3
   libNUMA enabled : $have_numa
   Using tcmalloc  : $have_tcmalloc
   Using jemalloc  : $have_jemalloc
   CPU profiler    : $have_profiler

   Hydro scheme       : $with_hydro
   Dimensionality     : $with_dimension
   Kernel function    : $with_kernel
   Equation of state  : $with_eos
   Adiabatic index    : $with_gamma
   Riemann solver     : $with_riemann
   Cooling function   : $with_cooling
   
   External potential  : $with_potential
   Multipole order     : $with_multipole_order
   No gravity below ID : $no_gravity_below_id

   Individual timers    : $enable_timers
   Task debugging       : $enable_task_debugging
   Threadpool debugging : $enable_threadpool_debugging
   Debugging checks     : $enable_debugging_checks
   Gravity checks       : $gravity_force_checks

 ------------------------])
