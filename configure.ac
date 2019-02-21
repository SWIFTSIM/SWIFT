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
AC_INIT([SWIFT],[0.8.0],[https://gitlab.cosma.dur.ac.uk/swift/swiftsim])
swift_config_flags="$*"

#  We want to stop when given unrecognised options. No subdirs so this is safe.
enable_option_checking=${enable_option_checking:-fatal}
if test -n "$ac_unrecognized_opts"; then
    case $enable_option_checking in
        no)
        ;;
        fatal)
            { $as_echo "$as_me: error: unrecognized options: $ac_unrecognized_opts" >&2
              { (exit 1); exit 1; }; }
        ;;
        *)
            $as_echo "$as_me: WARNING: unrecognized options: $ac_unrecognized_opts" >&2
        ;;
    esac
fi

AC_COPYRIGHT
AC_CONFIG_SRCDIR([src/space.c])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([subdir-objects])

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

# If debug is selected then we also define SWIFT_DEVELOP_MODE to control
# any developer code options.
if test "x$ax_enable_debug" != "xno"; then
   AC_DEFINE([SWIFT_DEVELOP_MODE],1,[Enable developer code options])
fi

# Enable POSIX and platform extension preprocessor macros.
AC_USE_SYSTEM_EXTENSIONS

# Check for compiler version and vendor.
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

#  Restrict support.
AC_C_RESTRICT

# logger
AC_ARG_ENABLE([logger],
	[AS_HELP_STRING([--enable-logger],
		[enable the particle logger]
	)],
	[with_logger="${enableval}"],
	[with_logger="no"]
)

if test "$with_logger" = "yes"; then
   AC_DEFINE([WITH_LOGGER], 1, [logger enabled])
fi

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
      AX_COMPARE_VERSION($ax_cv_c_compiler_version, [ge], [5.0.0],
                          [
      : ${AR="gcc-ar"}
      : ${RANLIB="gcc-ranlib"}
                          ], [:] )
      AC_MSG_RESULT([added GCC interprocedural optimization support])
   elif test "$ax_cv_c_compiler_vendor" = "clang"; then
      CFLAGS="$CFLAGS -flto -fuse-ld=gold"
      LDFLAGS="$LDFLAGS -XCClinker -fuse-ld=gold"
      : ${AR="llvm-ar"}
      : ${LD="ld.gold"}
      : ${RANLIB="llvm-ranlib"}
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

# Indicate that MPIRUN can be modified by an environment variable
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
     [Store extra information for generating task dump files @<:@yes/no@:>@]
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

# Check if using our custom icbrtf is enalbled.
AC_ARG_ENABLE([custom-icbrtf],
   [AS_HELP_STRING([--enable-custom-icbrtf],
     [Use SWIFT's custom icbrtf function instead of the system cbrtf @<:@yes/no@:>@]
   )],
   [enable_custom_icbrtf="$enableval"],
   [enable_custom_icbrtf="no"]
)
if test "$enable_custom_icbrtf" = "yes"; then
   AC_DEFINE([WITH_ICBRTF],1,[Enable custom icbrtf])
fi

# Check whether we want to default to naive cell interactions
AC_ARG_ENABLE([naive-interactions],
   [AS_HELP_STRING([--enable-naive-interactions],
     [Activate use of naive cell interaction functions @<:@yes/no@:>@]
   )],
   [enable_naive_interactions="$enableval"],
   [enable_naive_interactions="no"]
)
if test "$enable_naive_interactions" = "yes"; then
   AC_DEFINE([SWIFT_USE_NAIVE_INTERACTIONS],1,[Enable use of naive cell interaction functions])
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

# Check whether we want to switch on glass making
AC_ARG_ENABLE([glass-making],
   [AS_HELP_STRING([--enable-glass-making],
     [Activate the glass-making procedure by reversing the sign of gravity @<:@yes/no@:>@]
   )],
   [gravity_glass_making="$enableval"],
   [gravity_glass_making="no"]
)
if test "$gravity_glass_making" == "yes"; then
   AC_DEFINE([SWIFT_MAKE_GRAVITY_GLASS], 1, [Make the code run in a way to produce a glass file for gravity/cosmology])
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

# Check whether we have any of the ARM v8.1 tick timers
AX_ASM_ARM_PMCCNTR
AX_ASM_ARM_CNTVCT

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

#  Disable hand written vectorisation. Slightly odd implementation as want
# to describe as --disable-hand-vec, but macro is enable (there is no enable action).
AC_ARG_ENABLE([hand-vec],
   [AS_HELP_STRING([--disable-hand-vec],
     [Disable intrinsic vectorization]
   )],
   [enable_hand_vec="$enableval"],
   [enable_hand_vec="yes"]
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
   elif test "$enable_hand_vec" = "yes"; then
      AC_DEFINE([WITH_VECTORIZATION],1,[Enable hand-written vectorization])
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

# Add the undefined sanitizer option to flags. Only useful for GCC
# version 4.9 and later and clang to detected undefined code behaviour
# such as integer overflow and memory alignment issues.
AC_ARG_ENABLE([undefined-sanitizer],
   [AS_HELP_STRING([--enable-undefined-sanitizer],
     [Enable detection of code that causes undefined behaviour @<:@no/yes@:>@]
   )],
   [enable_ubsan="$enableval"],
   [enable_ubsan="no"]
)

if test "$enable_ubsan" = "yes"; then
   if test "$ax_cv_c_compiler_vendor" = "gnu"; then
      AX_COMPARE_VERSION( $ax_cv_c_compiler_version, [ge], [4.9.0],
                          [enable_ubsan="yes"], [enable_ubsan="no"] )
   elif test "$ax_cv_c_compiler_vendor" = "clang"; then
      AX_COMPARE_VERSION( $ax_cv_c_compiler_version, [ge], [3.7.0],
                          [enable_ubsan="yes"], [enable_ubsan="no"] )
   fi
   if test "$enable_ubsan" = "yes"; then
      CFLAGS="$CFLAGS -fsanitize=undefined"
      AC_MSG_RESULT([added undefined sanitizer support])
   else
      AC_MSG_WARN([Compiler does not support undefined sanitizer option])
   fi
fi

# Autoconf stuff.
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_HEADER_STDC

# Check for the libraries we will need.
AC_CHECK_LIB(m,sqrt,,AC_MSG_ERROR(something is wrong with the math library!))

# Check for GSL. We test for this in the standard directories by default,
# and only disable if using --with-gsl=no or --without-gsl. When a value
# is given GSL must be found.
have_gsl="no"
AC_ARG_WITH([gsl],
    [AS_HELP_STRING([--with-gsl=PATH],
       [root directory where GSL is installed @<:@yes/no@:>@]
    )],
    [with_gsl="$withval"],
    [with_gsl="test"]
)
if test "x$with_gsl" != "xno"; then
   if test "x$with_gsl" != "xyes" -a "x$with_gsl" != "xtest" -a "x$with_gsl" != "x"; then
      GSL_LIBS="-L$with_gsl/lib -lgsl -lgslcblas"
      GSL_INCS="-I$with_gsl/include"
   else
      GSL_LIBS="-lgsl -lgslcblas"
      GSL_INCS=""
   fi
   #  GSL is not specified, so just check if we have it.
   if test "x$with_gsl" = "xtest"; then
      AC_CHECK_LIB([gslcblas],[cblas_dgemm],[have_gsl="yes"],[have_gsl="no"],$GSL_LIBS)
      if test "x$have_gsl" != "xno"; then
         AC_DEFINE([HAVE_LIBGSLCBLAS],1,[The GSL CBLAS library appears to be present.])
         AC_CHECK_LIB([gsl],[gsl_integration_qag],
            AC_DEFINE([HAVE_LIBGSL],1,[The GSL library appears to be present.]),
            [have_gsl="no"],$GSL_LIBS)
      fi
   else
      AC_CHECK_LIB([gslcblas],[cblas_dgemm],
         AC_DEFINE([HAVE_LIBGSLCBLAS],1,[The GSL CBLAS library appears to be present.]),
         AC_MSG_ERROR(something is wrong with the GSL CBLAS library!), $GSL_LIBS)
      AC_CHECK_LIB([gsl],[gsl_integration_qag],
         AC_DEFINE([HAVE_LIBGSL],1,[The GSL library appears to be present.]),
         AC_MSG_ERROR(something is wrong with the GSL library!), $GSL_LIBS)
      have_gsl="yes"
   fi
   if test "$have_gsl" = "no"; then
      GSL_LIBS=""
      GSL_INCS=""
   fi
fi
AC_SUBST([GSL_LIBS])
AC_SUBST([GSL_INCS])
AM_CONDITIONAL([HAVEGSL],[test -n "$GSL_LIBS"])

# Check for pthreads.
AX_PTHREAD([LIBS="$PTHREAD_LIBS $LIBS" CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
    CC="$PTHREAD_CC" LDFLAGS="$LDFLAGS $PTHREAD_LIBS $LIBS"],
    AC_MSG_ERROR([Could not find a working version of
    the pthread library. Make sure you have the library and header files installed
    or use CPPFLAGS and LDFLAGS if the library is installed in a
    non-standard location.]))

# Check whether POSIX thread barriers are implemented (e.g. OSX does not have them)
have_pthread_barrier="no"
AC_CHECK_LIB(pthread, pthread_barrier_init,
	     have_pthread_barrier="yes",
	     AC_MSG_WARN(POSIX implementation does not have barriers. SWIFT will use home-made ones.))
if test "x$have_pthread_barrier" == "xyes"; then
  AC_DEFINE([HAVE_PTHREAD_BARRIERS], [1], [The posix library implements barriers])
fi

# Check whether POSIX file allocation functions exist (e.g. OSX does not have them)
AC_CHECK_LIB(pthread, posix_fallocate,
	     AC_DEFINE([HAVE_POSIX_FALLOCATE], [1], [The posix library implements file allocation functions.]),
	     AC_MSG_WARN(POSIX implementation does not have file allocation functions.))

# Check for METIS.
have_metis="no"
AC_ARG_WITH([metis],
    [AS_HELP_STRING([--with-metis=PATH],
       [root directory where METIS is installed @<:@yes/no@:>@]
    )],
    [with_metis="$withval"],
    [with_metis="no"]
)

METIS_LIBS=""
if test "x$with_metis" != "xno"; then

# Check if we have METIS.
   if test "x$with_metis" != "xyes" -a "x$with_metis" != "x"; then
      METIS_LIBS="-L$with_metis/lib -lmetis"
      METIS_INCS="-I$with_metis/include"
   else
      METIS_LIBS="-lmetis"
      METIS_INCS=""
   fi
   AC_CHECK_LIB([metis],[METIS_PartGraphKway], [have_metis="yes"],
                [have_metis="no"], $METIS_LIBS)
   if test "$have_metis" == "yes"; then
      AC_DEFINE([HAVE_METIS],1,[The METIS library is present.])
   else
      AC_MSG_ERROR("Failed to find a METIS library")
   fi
fi

AC_SUBST([METIS_LIBS])
AC_SUBST([METIS_INCS])
AM_CONDITIONAL([HAVEMETIS],[test -n "$METIS_LIBS"])

# Check for ParMETIS note we can have both as ParMETIS uses METIS.
have_parmetis="no"
AC_ARG_WITH([parmetis],
    [AS_HELP_STRING([--with-parmetis=PATH],
       [root directory where ParMETIS is installed @<:@yes/no@:>@]
    )],
    [with_parmetis="$withval"],
    [with_parmetis="no"]
)

if test "x$with_parmetis" != "xno"; then

# Check if we have ParMETIS.
   if test "x$with_parmetis" != "xyes" -a "x$with_parmetis" != "x"; then
      PARMETIS_LIBS="-L$with_parmetis/lib -lparmetis"
      PARMETIS_INCS="-I$with_parmetis/include"
   else
      PARMETIS_LIBS="-lparmetis"
      PARMETIS_INCS=""
   fi
   AC_CHECK_LIB([parmetis],[ParMETIS_V3_RefineKway], [have_parmetis="yes"],
                [have_parmetis="no"], $PARMETIS_LIBS)
   if test "$have_parmetis" == "no"; then

# A build may use an external METIS library, check for that.

      if test "x$with_parmetis" != "xyes" -a "x$with_parmetis" != "x"; then
         PARMETIS_LIBS="-L$with_parmetis/lib -lparmetis -lmetis"
         PARMETIS_INCS="-I$with_parmetis/include"
      else
         PARMETIS_LIBS="-lparmetis -lmetis"
         PARMETIS_INCS=""
      fi
      AC_CHECK_LIB([parmetis],[ParMETIS_V3_RefineKway], [have_parmetis="yes"],
                   [have_parmetis="no"], [$METIS_LIBS $PARMETIS_LIBS])

   fi
   if test "$have_parmetis" == "yes"; then
      AC_DEFINE([HAVE_PARMETIS],1,[The ParMETIS library is present.])
   else
      AC_MSG_ERROR("Failed to find a ParMETIS library")
   fi
fi

AC_SUBST([PARMETIS_LIBS])
AC_SUBST([PARMETIS_INCS])
AM_CONDITIONAL([HAVEPARMETIS],[test -n "$PARMETIS_LIBS"])

# METIS fixed width integer printing can require this, so define. Only needed
# for some non C99 compilers, i.e. C++ pre C++11.
AH_VERBATIM([__STDC_FORMAT_MACROS],
            [/* Needed to get PRIxxx macros from stdint.h when not using C99 */
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS 1
#endif])

# Check for FFTW. We test for this in the standard directories by default,
# and only disable if using --with-fftw=no or --without-fftw. When a value
# is given FFTW must be found.
# If FFTW is found, we check whether this is the threaded version.
have_fftw="no"
AC_ARG_WITH([fftw],
    [AS_HELP_STRING([--with-fftw=PATH],
       [root directory where fftw is installed @<:@yes/no@:>@]
    )],
    [with_fftw="$withval"],
    [with_fftw="test"]
)
if test "x$with_fftw" != "xno"; then

   # Was FFTW's location specifically given?
   if test "x$with_fftw" != "xyes" -a "x$with_fftw" != "xtest" -a "x$with_fftw" != "x"; then
      FFTW_LIBS="-L$with_fftw/lib -lfftw3"
      FFTW_INCS="-I$with_fftw/include"
   else
      FFTW_LIBS="-lfftw3"
      FFTW_INCS=""
   fi

   #  FFTW is not specified, so just check if we have it.
   if test "x$with_fftw" = "xtest"; then
      AC_CHECK_LIB([fftw3],[fftw_malloc],[have_fftw="yes"],[have_fftw="no"],$FFTW_LIBS)
      if test "x$have_fftw" != "xno"; then
      	 AC_DEFINE([HAVE_FFTW],1,[The FFTW library appears to be present.])
      fi
   # FFTW was specified, check that it was a valid location.
   else
      AC_CHECK_LIB([fftw3],[fftw_malloc],
         AC_DEFINE([HAVE_FFTW],1,[The FFTW library appears to be present.]),
         AC_MSG_ERROR(something is wrong with the FFTW library!), $FFTW_LIBS)
      have_fftw="yes"
   fi

   # FFTW was requested not to be used.
   if test "$have_fftw" = "no"; then
      FFTW_LIBS=""
      FFTW_INCS=""
   fi

   # Now, check whether we have the threaded version of FFTW
   if test "x$have_fftw" = "xyes"; then

      # Was FFTW's location specifically given?
      if test "x$with_fftw" != "xyes" -a "x$with_fftw" != "xtest" -a "x$with_fftw" != "x"; then
        FFTW_THREADED_LIBS="-L$with_fftw/lib -lfftw3_threads -lfftw3"
        FFTW_THREADED_INCS="-I$with_fftw/include"
      else
        FFTW_THREADED_LIBS="-lfftw3_threads -lfftw3"
        FFTW_THREADED_INCS=""
      fi

      # Verify that the library is threaded
      AC_CHECK_LIB([fftw3],[fftw_init_threads],[have_threaded_fftw="yes"],
		   [have_threaded_fftw="no"], $FFTW_THREADED_LIBS)

      # If found, update things
      if test "x$have_threaded_fftw" = "xyes"; then
         AC_DEFINE([HAVE_THREADED_FFTW],1,[The threaded FFTW library appears to be present.])
         FFTW_LIBS=$FFTW_THREADED_LIBS
         FFTW_INCS=$FFTW_THREADED_INCS
	 have_fftw="yes - threaded"
      fi
   fi
fi

AC_ARG_WITH([arm-fftw],
    [AS_HELP_STRING([--with-arm-fftw=PATH],
      [root directory where arm fft library is installed @<:@yes/no@:>@]
    )],
    [with_arm_fftw="$withval"],
    [with_arm_fftw=no]
)
if test "x$with_arm_fftw" != "xno"; then

   # Was FFTW's location specifically given?
   if test "x$with_arm_fftw" != "xyes" -a "x$with_arm_fftw" != "xtest" -a "x$with_arm_fftw" != "x"; then
      FFTW_LIBS="-L$with_arm_fftw/lib -larmpl_lp64"
      FFTW_INCS="-I$with_arm_fftw/include"
   else
      FFTW_LIBS="-larmpl_lp64"
      FFTW_INCS=""
   fi

   #  FFTW is not specified, so just check if we have it.
   if test "x$with_arm_fftw" = "xtest"; then
      AC_CHECK_LIB([armpl_lp64],[fftw_malloc],[have_fftw="yes"],[have_fftw="no"],$FFTW_LIBS)
      if test "x$have_arm_fftw" != "xno"; then
      	 AC_DEFINE([HAVE_FFTW],1,[The FFTW library appears to be present.])
	 have_fftw="yes - ARM"
      fi
   # FFTW was specified, check that it was a valid location.
   else
      AC_CHECK_LIB([armpl_lp64],[fftw_malloc],
         AC_DEFINE([HAVE_FFTW],1,[The FFTW library appears to be present.]),
         AC_MSG_ERROR(something is wrong with the FFTW library!), $FFTW_LIBS)
      have_fftw="yes - ARM"
   fi

   # FFTW was requested not to be used.
   if test "$have_arm_fftw" = "no"; then
      FFTW_LIBS=""
      FFTW_INCS=""
   fi

   # Now, check whether we have the threaded version of FFTW
   if test "x$have_arm_fftw" = "xyes"; then

      # Was FFTW's location specifically given?
      if test "x$with_arm_fftw" != "xyes" -a "x$with_arm_fftw" != "xtest" -a "x$with_arm_fftw" != "x"; then
        FFTW_THREADED_LIBS="-L$with_arm_fftw/lib -larmpl_lp64_threads -larmpl_lp64"
        FFTW_THREADED_INCS="-I$with_arm_fftw/include"
      else
        FFTW_THREADED_LIBS="-larmpl_lp64_threads -larmpl_lp64"
        FFTW_THREADED_INCS=""
      fi

      # Verify that the library is threaded
      AC_CHECK_LIB([armpl_lp64],[fftw_init_threads],[have_threaded_fftw="yes"],
                  [have_threaded_fftw="no"], $FFTW_THREADED_LIBS)

      # If found, update things
      if test "x$have_threaded_fftw" = "xyes"; then
         AC_DEFINE([HAVE_THREADED_FFTW],1,[The threaded FFTW library appears to be present.])
         FFTW_LIBS=$FFTW_THREADED_LIBS
         FFTW_INCS=$FFTW_THREADED_INCS
         have_fftw="yes - ARM - threaded"
      fi
   fi
fi
AC_SUBST([FFTW_LIBS])
AC_SUBST([FFTW_INCS])
AM_CONDITIONAL([HAVEFFTW],[test -n "$FFTW_LIBS"])

#  Check for -lprofiler usually part of the gperftools along with tcmalloc.
have_profiler="no"
AC_ARG_WITH([profiler],
   [AS_HELP_STRING([--with-profiler=PATH],
      [use cpu profiler library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_profiler="$withval"],
   [with_profiler="no"]
)
if test "x$with_profiler" != "xno"; then
   if test "x$with_profiler" != "xyes" -a "x$with_profiler" != "x"; then
      proflibs="-L$with_profiler -lprofiler"
   else
      proflibs="-lprofiler"
   fi
   AC_CHECK_LIB([profiler],[ProfilerFlush],
    [have_profiler="yes"
      AC_DEFINE([WITH_PROFILER],1,[Link against the gperftools profiling library.])],
    [have_profiler="no"], $proflibs)

   if test "$have_profiler" = "yes"; then
      PROFILER_LIBS="$proflibs"
   else
      PROFILER_LIBS=""
   fi
fi
AC_SUBST([PROFILER_LIBS])
AM_CONDITIONAL([HAVEPROFILER],[test -n "$PROFILER_LIBS"])

# Check for special allocators
have_special_allocator="no"

#  Check for tcmalloc a fast malloc that is part of the gperftools.
have_tcmalloc="no"
AC_ARG_WITH([tcmalloc],
   [AS_HELP_STRING([--with-tcmalloc=PATH],
      [use tcmalloc library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_tcmalloc="$withval"],
   [with_tcmalloc="no"]
)
if test "x$with_tcmalloc" != "xno" -a "x$have_special_allocator" != "xno"; then
   AC_MSG_ERROR("Cannot activate more than one alternative malloc library")
fi

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

      AC_DEFINE([HAVE_TCMALLOC],1,[The tcmalloc library appears to be present.])

      have_special_allocator="tcmalloc"

      # Prevent compilers that replace the calls with built-ins (GNU 99) from doing so.
      case "$ax_cv_c_compiler_vendor" in
        intel | gnu | clang)
             CFLAGS="$CFLAGS -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free"
          ;;
      esac

   else
      TCMALLOC_LIBS=""
   fi
fi
AC_SUBST([TCMALLOC_LIBS])
AM_CONDITIONAL([HAVETCMALLOC],[test -n "$TCMALLOC_LIBS"])

#  Check for jemalloc another fast malloc that is good with contention.
have_jemalloc="no"
AC_ARG_WITH([jemalloc],
   [AS_HELP_STRING([--with-jemalloc=PATH],
      [use jemalloc library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_jemalloc="$withval"],
   [with_jemalloc="no"]
)
if test "x$with_jemalloc" != "xno" -a "x$have_special_allocator" != "xno"; then
   AC_MSG_ERROR("Cannot activate more than one alternative malloc library")
fi

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

      AC_DEFINE([HAVE_JEMALLOC],1,[The jemalloc library appears to be present.])

      have_special_allocator="jemalloc"

      # Prevent compilers that replace the regular calls with built-ins (GNU 99) from doing so.
      case "$ax_cv_c_compiler_vendor" in
        intel | gnu | clang)
             CFLAGS="$CFLAGS -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free"
          ;;
      esac

   else
      JEMALLOC_LIBS=""
   fi
fi
AC_SUBST([JEMALLOC_LIBS])
AM_CONDITIONAL([HAVEJEMALLOC],[test -n "$JEMALLOC_LIBS"])

#  Check for tbbmalloc, Intel's fast and parallel allocator
have_tbbmalloc="no"
AC_ARG_WITH([tbbmalloc],
   [AS_HELP_STRING([--with-tbbmalloc=PATH],
      [use tbbmalloc library or specify the directory with lib @<:@yes/no@:>@]
   )],
   [with_tbbmalloc="$withval"],
   [with_tbbmalloc="no"]
)
if test "x$with_tbbmalloc" != "xno" -a "x$have_special_allocator" != "xno"; then
   AC_MSG_ERROR("Cannot activate more than one alternative malloc library")
fi

if test "x$with_tbbmalloc" != "xno"; then
   if test "x$with_tbbmalloc" != "xyes" -a "x$with_tbbmalloc" != "x"; then
      tbblibs="-L$with_tbbmalloc -ltbbmalloc_proxy -ltbbmalloc"
   else
      tbblibs="-ltbbmalloc_proxy -ltbbmalloc"
   fi
   AC_CHECK_LIB([tbbmalloc],[scalable_malloc],[have_tbbmalloc="yes"],[have_tbbmalloc="no"],
                $tbblibs)

   if test "$have_tbbmalloc" = "yes"; then
      TBBMALLOC_LIBS="$tbblibs"

      AC_DEFINE([HAVE_TBBMALLOC],1,[The TBBmalloc library appears to be present.])

      have_special_allocator="TBBmalloc"

      # Prevent compilers that replace the calls with built-ins (GNU 99) from doing so.
      case "$ax_cv_c_compiler_vendor" in
        intel | gnu | clang)
             CFLAGS="$CFLAGS -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free"
          ;;
      esac

   else
      TBBMALLOC_LIBS=""
   fi
fi
AC_SUBST([TBBMALLOC_LIBS])
AM_CONDITIONAL([HAVETBBMALLOC],[test -n "$TBBMALLOC_LIBS"])

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
        old_CPPFLAGS="$CPPFLAGS"
        CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
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
        CPPFLAGS="$old_CPPFLAGS"
    fi
fi
AM_CONDITIONAL([HAVEPARALLELHDF5],[test "$have_parallel_hdf5" = "yes"])

# Check for grackle.
have_grackle="no"
AC_ARG_WITH([grackle],
    [AS_HELP_STRING([--with-grackle=PATH],
       [root directory where grackle is installed @<:@yes/no@:>@]
    )],
    [with_grackle="$withval"],
    [with_grackle="no"]
)
if test "x$with_grackle" != "xno"; then
   AC_PROG_FC
   AC_FC_LIBRARY_LDFLAGS
   if test "x$with_grackle" != "xyes" -a "x$with_grackle" != "x"; then
      GRACKLE_LIBS="-L$with_grackle/lib -lgrackle"
      GRACKLE_INCS="-I$with_grackle/include"
   else
      GRACKLE_LIBS="-lgrackle"
      GRACKLE_INCS=""
   fi

   have_grackle="yes"

   echo $GRACKLE_LIBS

   AC_CHECK_LIB(
      [grackle],
      [initialize_chemistry_data],
      [AC_DEFINE([HAVE_GRACKLE],1,[The GRACKLE library appears to be present.])
        AC_DEFINE([CONFIG_BFLOAT_8],1,[Use doubles in grackle])
      ],
      [AC_MSG_ERROR(Cannot find grackle library!)],
      [$GRACKLE_LIBS])
fi
AC_SUBST([GRACKLE_LIBS])
AC_SUBST([GRACKLE_INCS])
AM_CONDITIONAL([HAVEGRACKLE],[test -n "$GRACKLE_LIBS"])

# Check for VELOCIraptor.
have_velociraptor="no"
AC_ARG_WITH([velociraptor],
    [AS_HELP_STRING([--with-velociraptor=PATH],
       [Directory where velociraptor library exists @<:@yes/no@:>@]
    )],
    [with_velociraptor="$withval"],
    [with_velociraptor="no"]
)
if test "x$with_velociraptor" != "xno"; then
   AC_PROG_FC
   AC_FC_LIBRARY_LDFLAGS
   if test "x$with_velociraptor" != "xyes" -a "x$with_velociraptor" != "x"; then
      VELOCIRAPTOR_LIBS="-L$with_velociraptor -lvelociraptor -lmpi -lstdc++ -lhdf5_cpp"
      CFLAGS="$CFLAGS -fopenmp"
   else
      VELOCIRAPTOR_LIBS=""
   fi

   have_velociraptor="yes"

   AC_CHECK_LIB(
      [velociraptor],
      [InitVelociraptor],
      [AC_DEFINE([HAVE_VELOCIRAPTOR],1,[The VELOCIraptor library appears to be present.])],
      [AC_MSG_ERROR(Cannot find VELOCIraptor library at $with_velociraptor or incompatible HDF5 library loaded.)],
      [$VELOCIRAPTOR_LIBS $HDF5_LDFLAGS $HDF5_LIBS $GSL_LIBS]
   )
fi
AC_SUBST([VELOCIRAPTOR_LIBS])
AM_CONDITIONAL([HAVEVELOCIRAPTOR],[test -n "$VELOCIRAPTOR_LIBS"])

# Check for dummy VELOCIraptor.
AC_ARG_ENABLE([dummy-velociraptor],
    [AS_HELP_STRING([--enable-dummy-velociraptor],
       [Enable dummy velociraptor compilation @<:@yes/no@:>@]
    )],
    [enable_dummy_velociraptor="$enableval"],
    [enable_dummy_velociraptor="no"]
)

if test "$enable_dummy_velociraptor" = "yes"; then
  have_velociraptor="yes"

  AC_DEFINE(HAVE_VELOCIRAPTOR,1,[The VELOCIraptor library appears to be present.])
  AC_DEFINE(HAVE_DUMMY_VELOCIRAPTOR,1,[The dummy VELOCIraptor library is present.])
fi

# Check for floating-point execeptions
AC_CHECK_FUNC(feenableexcept, AC_DEFINE([HAVE_FE_ENABLE_EXCEPT],[1],
    [Defined if the floating-point exception can be enabled using non-standard GNU functions.]))

# Check for setaffinity.
AC_CHECK_FUNC(pthread_setaffinity_np, AC_DEFINE([HAVE_SETAFFINITY],[1],
    [Defined if pthread_setaffinity_np exists.]) )
AM_CONDITIONAL(HAVESETAFFINITY,
    [test "$ac_cv_func_pthread_setaffinity_np" = "yes"])

# If available check for NUMA as well. There is a problem with the headers of
# this library, mainly that they do not pass the strict prototypes check when
# installed outside of the system directories. So we actually do this check
# in two phases. The basic ones first (before strict-prototypes is added to CFLAGS).
have_numa="no"
AC_ARG_WITH([numa],
    [AS_HELP_STRING([--with-numa=PATH],
       [Directory where the NUMA library exists @<:@yes/no@:>@]
    )],
    [with_numa="$withval"],
    [with_numa="yes"]
)
if test "$ac_cv_func_pthread_setaffinity_np" = "yes" -a "x$with_numa" != "xno"; then

    if test "x$with_numa" != "xyes" -a "x$with_numa" != "x"; then
        NUMA_LIBS="-L$with_numa/lib -lnuma"
        NUMA_INCS="-I$with_numa/include"
    else
        NUMA_LIBS="-lnuma"
        NUMA_INCS=""
    fi

    #  Test for header file.
    old_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $NUMA_INCS"
    AC_CHECK_HEADER([numa.h])
    CPPFLAGS="$old_CPPFLAGS"
    if test "$ac_cv_header_numa_h" = "yes"; then

        #  If NUMA location is specified check if we have it.
        if test "x$with_numa" != "xyes" -a "x$with_numa" != "x"; then
            AC_CHECK_LIB([numa],[numa_available],
                AC_DEFINE([HAVE_LIBNUMA],1,[The NUMA library appears to be present.]),
                AC_MSG_ERROR(something is wrong with the NUMA library!), $NUMA_LIBS)
            have_numa="yes"
        else
            AC_CHECK_LIB([numa],[numa_available],[have_numa="yes"],[have_numa="no"],$NUMA_LIBS)
            if test "x$have_numa" != "xno"; then
                AC_DEFINE([HAVE_LIBNUMA],1,[The NUMA library appears to be present.])
            fi
        fi
    fi

    #  We can live without this.
    if test "$have_numa" = "no"; then
       NUMA_LIBS=""
    fi
fi
AC_SUBST([NUMA_LIBS])

# Check for Intel and PowerPC intrinsics header optionally used by vector.h.
AC_CHECK_HEADERS([immintrin.h], [], [],
[#ifdef HAVE_IMMINTRIN_H
# include <immintrin.h>
#endif
])
AC_CHECK_HEADERS([altivec.h], [], [],
[#ifdef HAVE_ALTIVEC_H
# include <altivec.h>
#endif
])

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

# Special timers for the ARM v7 platforms (taken from FFTW-3 to match their cycle.h)
AC_ARG_ENABLE(armv7a-cntvct, [AC_HELP_STRING([--enable-armv7a-cntvct],[enable the cycle counter on Armv7a via the CNTVCT register])], have_armv7acntvct=$enableval)
if test "$have_armv7acntvct"x = "yes"x; then
	AC_DEFINE(HAVE_ARMV7A_CNTVCT,1,[Define if you have enabled the CNTVCT cycle counter on ARMv7a])
fi

AC_ARG_ENABLE(armv7a-pmccntr, [AC_HELP_STRING([--enable-armv7a-pmccntr],[enable the cycle counter on Armv7a via the PMCCNTR register])], have_armv7apmccntr=$enableval)
if test "$have_armv7apmccntr"x = "yes"x; then
	AC_DEFINE(HAVE_ARMV7A_PMCCNTR,1,[Define if you have enabled the PMCCNTR cycle counter on ARMv7a])
fi

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
             CFLAGS="$CFLAGS -Wall -Wextra -Wno-unused-parameter -Wshadow"
          ;;
	  intel)
             CFLAGS="$CFLAGS -w2 -Wunused-variable -Wshadow"
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

    # We want strict-prototypes, but this must still work even if warnings
    # are an error.
    AX_CHECK_COMPILE_FLAG([-Wstrict-prototypes],[CFLAGS="$CFLAGS -Wstrict-prototypes"],
                          [CFLAGS="$CFLAGS"],[$CFLAGS],[AC_LANG_SOURCE([int main(void){return 0;}])])
fi

# Second part of the NUMA library checks. We now decide if we need to use
# -isystem to get around the strict-prototypes problem. Assumes isystem
# is available when strict-prototypes is.
if test "$have_numa" != "no"; then
    if test "x$with_numa" != "xyes" -a "x$with_numa" != "x"; then
        case "$CFLAGS" in
            *strict-prototypes*)
                NUMA_INCS="-isystem$with_numa/include"
                # This may still fail if CPATH is used, so we check if the
                # headers are usable.
                AS_UNSET(ac_cv_header_numa_h)
                old_CPPFLAGS="$CPPFLAGS"
                CPPFLAGS="$CPPFLAGS $NUMA_INCS"
                numa_failed="no"
                AC_CHECK_HEADER([numa.h],[numa_failed="no"],
                                [numa_failed="yes"])
                if test "$numa_failed" = "yes"; then
                    AC_MSG_ERROR([Failed to compile the numa.h header file: you may need to set --enable-compiler-warnings to yes or no])
                fi
                CPPFLAGS="$old_CPPFLAGS"
            ;;
            *)
                NUMA_INCS="-I$with_numa/include"
            ;;
        esac
   fi
fi
AC_SUBST([NUMA_INCS])

# Various package configuration options.

# Master subgrid options
# If you add a restriction (e.g. no cooling, chemistry or hydro)
# you will need to check for overwrite after reading the additional options.
# As an example for this, see the call to AC_ARG_WITH for cooling.
AC_ARG_WITH([subgrid],
	[AS_HELP_STRING([--with-subgrid=<subgrid>],
		[Master switch for subgrid methods. Inexperienced user should start from here @<:@none, GEAR, EAGLE default: none@:>@]
	)],
	[with_subgrid="$withval"],
	[with_subgrid=none]
)

# Default values
with_subgrid_cooling=none
with_subgrid_chemistry=none
with_subgrid_tracers=none
with_subgrid_entropy_floor=none
with_subgrid_stars=none
with_subgrid_star_formation=none
with_subgrid_feedback=none

case "$with_subgrid" in
   yes)
      AC_MSG_ERROR([Invalid option. A subgrid model must be chosen.])
   ;;
   none)
   ;;
   GEAR)
	with_subgrid_cooling=grackle
	with_subgrid_chemistry=GEAR
	with_subgrid_tracers=none
	with_subgrid_entropy_floor=none
	with_subgrid_stars=GEAR
	with_subgrid_star_formation=none
	with_subgrid_feedback=thermal
   ;;
   EAGLE)
	with_subgrid_cooling=EAGLE
	with_subgrid_chemistry=EAGLE
	with_subgrid_tracers=EAGLE
	with_subgrid_entropy_floor=EAGLE
	with_subgrid_stars=EAGLE
	with_subgrid_star_formation=EAGLE
	with_subgrid_feedback=none
   ;;
   *)
      AC_MSG_ERROR([Unknown subgrid choice: $with_subgrid])
   ;;
esac

# Gravity scheme.
AC_ARG_WITH([gravity],
   [AS_HELP_STRING([--with-gravity=<scheme>],
      [Gravity scheme to use @<:@default, with-potential, default: default@:>@]
   )],
   [with_gravity="$withval"],
   [with_gravity="default"]
)

case "$with_gravity" in
   with-potential)
      AC_DEFINE([POTENTIAL_GRAVITY], [1], [Gravity scheme with potential calculation])
   ;;
   default)
      AC_DEFINE([DEFAULT_GRAVITY], [1], [Default gravity scheme])
   ;;
   *)
      AC_MSG_ERROR([Unknown gravity scheme: $with_gravity])
   ;;
esac

# Hydro scheme.
AC_ARG_WITH([hydro],
   [AS_HELP_STRING([--with-hydro=<scheme>],
      [Hydro dynamics to use @<:@gadget2, minimal, pressure-entropy, pressure-energy, pressure-energy-monaghan, default, gizmo-mfv, gizmo-mfm, shadowfax, planetary, anarchy-pu debug default: gadget2@:>@]
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
   pressure-entropy)
      AC_DEFINE([HOPKINS_PE_SPH], [1], [Pressure-Entropy SPH])
   ;;
   pressure-energy)
      AC_DEFINE([HOPKINS_PU_SPH], [1], [Pressure-Energy SPH])
   ;;
   pressure-energy-monaghan)
      AC_DEFINE([HOPKINS_PU_SPH_MONAGHAN], [1], [Pressure-Energy SPH with M&M Variable A.V.])
   ;;
   default)
      AC_DEFINE([DEFAULT_SPH], [1], [Default SPH])
   ;;
   gizmo-mfv)
      AC_DEFINE([GIZMO_MFV_SPH], [1], [GIZMO MFV SPH])
   ;;
   gizmo-mfm)
      AC_DEFINE([GIZMO_MFM_SPH], [1], [GIZMO MFM SPH])
   ;;
   shadowfax)
      AC_DEFINE([SHADOWFAX_SPH], [1], [Shadowfax SPH])
   ;;
   planetary)
      AC_DEFINE([PLANETARY_SPH], [1], [Planetary SPH])
   ;;
   anarchy-pu)
      AC_DEFINE([ANARCHY_PU_SPH], [1], [ANARCHY (PU) SPH])
   ;;


   *)
      AC_MSG_ERROR([Unknown hydrodynamics scheme: $with_hydro])
   ;;
esac

# Check if debugging interactions stars is switched on.
AC_ARG_ENABLE([debug-interactions-stars],
   [AS_HELP_STRING([--enable-debug-interactions-stars],
     [Activate interaction debugging for stars, logging a maximum of @<:@N@:>@ neighbours. Defaults to 256 if no value set.]
   )],
   [enable_debug_interactions_stars="$enableval"],
   [enable_debug_interactions_stars="no"]
)
if test "$enable_debug_interactions_stars" != "no"; then
    AC_DEFINE([DEBUG_INTERACTIONS_STARS],1,[Enable interaction debugging for stars])
    if test "$enable_debug_interactions_stars" == "yes"; then
      AC_DEFINE([MAX_NUM_OF_NEIGHBOURS_STARS],256,[The maximum number of particle neighbours to be logged for stars])
      [enable_debug_interactions_stars="yes (Logging up to 256 neighbours)"]
    else
      AC_DEFINE_UNQUOTED([MAX_NUM_OF_NEIGHBOURS_STARS], [$enableval] ,[The maximum number of particle neighbours to be logged for stars])
      [enable_debug_interactions_stars="yes (Logging up to $enableval neighbours)"]
    fi
fi

# Check if debugging interactions is switched on.
AC_ARG_ENABLE([debug-interactions],
   [AS_HELP_STRING([--enable-debug-interactions],
     [Activate interaction debugging, logging a maximum of @<:@N@:>@ neighbours. Defaults to 256 if no value set.]
   )],
   [enable_debug_interactions="$enableval"],
   [enable_debug_interactions="no"]
)
if test "$enable_debug_interactions" != "no"; then
  if test "$with_hydro" = "gadget2"; then
      AC_DEFINE([DEBUG_INTERACTIONS_SPH],1,[Enable interaction debugging])
    if test "$enable_debug_interactions" == "yes"; then
      AC_DEFINE([MAX_NUM_OF_NEIGHBOURS],256,[The maximum number of particle neighbours to be logged])
      [enable_debug_interactions="yes (Logging up to 256 neighbours)"]
    else
      AC_DEFINE_UNQUOTED([MAX_NUM_OF_NEIGHBOURS], [$enableval] ,[The maximum number of particle neighbours to be logged])
      [enable_debug_interactions="yes (Logging up to $enableval neighbours)"]
    fi
  else
    [enable_debug_interactions="no (only available for gadget2 hydro scheme)"]
  fi
fi


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
      [equation of state @<:@ideal-gas, isothermal-gas, planetary default: ideal-gas@:>@]
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
   planetary)
      AC_DEFINE([EOS_PLANETARY], [1], [All planetary equations of state])
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
      [cooling function @<:@none, const-du, const-lambda, EAGLE, grackle, grackle1, grackle2, grackle3 default: none@:>@]
   )],
   [with_cooling="$withval"],
   [with_cooling="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_cooling" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-cooling together])
   else
      with_cooling="$with_subgrid_cooling"
   fi
fi

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
   compton)
      AC_DEFINE([COOLING_COMPTON], [1], [Compton cooling off the CMB])
   ;;
   grackle)
      AC_DEFINE([COOLING_GRACKLE], [1], [Cooling via the grackle library])
      AC_DEFINE([COOLING_GRACKLE_MODE], [0], [Grackle chemistry network, mode 0])
   ;;
   grackle1)
      AC_DEFINE([COOLING_GRACKLE], [1], [Cooling via the grackle library])
      AC_DEFINE([COOLING_GRACKLE_MODE], [1], [Grackle chemistry network, mode 1])
   ;;
   grackle2)
      AC_DEFINE([COOLING_GRACKLE], [1], [Cooling via the grackle library])
      AC_DEFINE([COOLING_GRACKLE_MODE], [2], [Grackle chemistry network, mode 2])
   ;;
   grackle3)
      AC_DEFINE([COOLING_GRACKLE], [1], [Cooling via the grackle library])
      AC_DEFINE([COOLING_GRACKLE_MODE], [3], [Grackle chemistry network, mode 3])
   ;;
   EAGLE)
      AC_DEFINE([COOLING_EAGLE], [1], [Cooling following the EAGLE model])
   ;;
   *)
      AC_MSG_ERROR([Unknown cooling function: $with_cooling])
   ;;
esac

#  chemistry function
AC_ARG_WITH([chemistry],
   [AS_HELP_STRING([--with-chemistry=<function>],
      [chemistry function @<:@none, GEAR, EAGLE default: none@:>@]
   )],
   [with_chemistry="$withval"],
   [with_chemistry="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_chemistry" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-chemistry together])
   else
      with_chemistry="$with_subgrid_chemistry"
   fi
fi

case "$with_chemistry" in
   none)
      AC_DEFINE([CHEMISTRY_NONE], [1], [No chemistry function])
   ;;
   GEAR)
      AC_DEFINE([CHEMISTRY_GEAR], [1], [Chemistry taken from the GEAR model])
   ;;
   EAGLE)
      AC_DEFINE([CHEMISTRY_EAGLE], [1], [Chemistry taken from the EAGLE model])
   ;;
   *)
      AC_MSG_ERROR([Unknown chemistry function: $with_chemistry])
   ;;
esac

#  Particle tracers
AC_ARG_WITH([tracers],
   [AS_HELP_STRING([--with-tracers=<function>],
      [chemistry function @<:@none, EAGLE default: none@:>@]
   )],
   [with_tracers="$withval"],
   [with_tracers="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_tracers" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-tracers together])
   else
      with_tracers="$with_subgrid_tracers"
   fi
fi

case "$with_tracers" in
   none)
      AC_DEFINE([TRACERS_NONE], [1], [No tracers function])
   ;;
   EAGLE)
      AC_DEFINE([TRACERS_EAGLE], [1], [Tracers taken from the EAGLE model])
   ;;
   *)
      AC_MSG_ERROR([Unknown tracers choice: $with_tracers])
   ;;
esac

# Stellar model.
AC_ARG_WITH([stars],
   [AS_HELP_STRING([--with-stars=<model>],
      [Stellar model to use @<:@none, EAGLE, GEAR, debug default: none@:>@]
   )],
   [with_stars="$withval"],
   [with_stars="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_stars" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-stars together])
   else
      with_stars="$with_subgrid_stars"
   fi
fi

case "$with_stars" in
   EAGLE)
      AC_DEFINE([STARS_EAGLE], [1], [EAGLE stellar model])
   ;;
   GEAR)
      AC_DEFINE([STARS_GEAR], [1], [GEAR stellar model])
   ;;
   none)
      AC_DEFINE([STARS_NONE], [1], [None stellar model])
   ;;

   *)
      AC_MSG_ERROR([Unknown stellar model: $with_stars])
   ;;
esac

# Feedback model
AC_ARG_WITH([feedback],
   [AS_HELP_STRING([--with-feedback=<model>],
      [Feedback model to use @<:@none, thermal, debug default: none@:>@]
   )],
   [with_feedback="$withval"],
   [with_feedback="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_feedback" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-feedback together])
   else
      with_feedback="$with_subgrid_feedback"
   fi
fi

case "$with_feedback" in
   thermal)
      AC_DEFINE([FEEDBACK_THERMAL], [1], [Thermal Blastwave])
   ;;
   none)
   ;;

   *)
      AC_MSG_ERROR([Unknown feedback model: $with_feedback])
   ;;
esac

#  External potential
AC_ARG_WITH([ext-potential],
   [AS_HELP_STRING([--with-ext-potential=<pot>],
      [external potential @<:@none, point-mass, point-mass-ring, point-mass-softened, isothermal, softened-isothermal, nfw, hernquist, disc-patch, sine-wave, default: none@:>@]
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
   hernquist)
      AC_DEFINE([EXTERNAL_POTENTIAL_HERNQUIST], [1], [Hernquist external potential])
   ;;
   nfw)
      AC_DEFINE([EXTERNAL_POTENTIAL_NFW], [1], [Navarro-Frenk-White external potential])
   ;;
   disc-patch)
      AC_DEFINE([EXTERNAL_POTENTIAL_DISC_PATCH], [1], [Disc-patch external potential])
   ;;
   sine-wave)
      AC_DEFINE([EXTERNAL_POTENTIAL_SINE_WAVE], [1], [Sine wave external potential in 1D])
   ;;
   point-mass-ring)
      AC_DEFINE([EXTERNAL_POTENTIAL_POINTMASS_RING], [1], [Point mass potential for Keplerian Ring (Hopkins 2015).])
   ;;
   point-mass-softened)
      AC_DEFINE([EXTERNAL_POTENTIAL_POINTMASS_SOFT], [1], [Softened point-mass potential with form 1/(r^2 + softening^2).])
   ;;
   *)
      AC_MSG_ERROR([Unknown external potential: $with_potential])
   ;;
esac

#  Entropy floor
AC_ARG_WITH([entropy-floor], 
    [AS_HELP_STRING([--with-entropy-floor=<floor>],
       [entropy floor @<:@none, EAGLE, default: none@:>@] 
    )],
    [with_entropy_floor="$withval"],
    [with_entropy_floor="none"]
)
if test "$with_subgrid" != "none"; then
   if test "$with_entropy_floor" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-entropy-floor together])
   else
      with_entropy_floor="$with_subgrid_entropy_floor"
   fi
fi

case "$with_entropy_floor" in
   none)
      AC_DEFINE([ENTROPY_FLOOR_NONE], [1], [No entropy floor])
   ;;
   EAGLE)
      AC_DEFINE([ENTROPY_FLOOR_EAGLE], [1], [EAGLE entropy floor])
   ;;
   *)
      AC_MSG_ERROR([Unknown entropy floor model])
   ;;
esac 

#  Star formation
AC_ARG_WITH([star-formation], 
    [AS_HELP_STRING([--with-star-formation=<sfm>],
       [star formation @<:@none, EAGLE, default: none@:>@] 
    )],
    [with_star_formation="$withval"],
    [with_star_formation="none"]
)
if test "$with_subgrid" != "none"; then
   if test "$with_star_formation" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-star-formation together])
   else
      with_star_formation="$with_subgrid_star_formation"
   fi
fi

case "$with_star_formation" in
   none)
      AC_DEFINE([STAR_FORMATION_NONE], [1], [No star formation])
   ;;
   EAGLE)
      AC_DEFINE([STAR_FORMATION_EAGLE], [1], [EAGLE star formation model (Schaye and Dalla Vecchia (2008))])
   ;;
   *)
      AC_MSG_ERROR([Unknown star formation model])
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

# Check if using EAGLE cooling
AM_CONDITIONAL([HAVEEAGLECOOLING], [test $with_cooling = "EAGLE"])

# Handle .in files.
AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile examples/Cooling/CoolingRates/Makefile doc/Makefile doc/Doxyfile tests/Makefile])
AC_CONFIG_FILES([argparse/Makefile tools/Makefile])
AC_CONFIG_FILES([tests/testReading.sh], [chmod +x tests/testReading.sh])
AC_CONFIG_FILES([tests/testActivePair.sh], [chmod +x tests/testActivePair.sh])
AC_CONFIG_FILES([tests/test27cells.sh], [chmod +x tests/test27cells.sh])
AC_CONFIG_FILES([tests/test27cellsPerturbed.sh], [chmod +x tests/test27cellsPerturbed.sh])
AC_CONFIG_FILES([tests/test27cellsStars.sh], [chmod +x tests/test27cellsStars.sh])
AC_CONFIG_FILES([tests/test27cellsStarsPerturbed.sh], [chmod +x tests/test27cellsStarsPerturbed.sh])
AC_CONFIG_FILES([tests/test125cells.sh], [chmod +x tests/test125cells.sh])
AC_CONFIG_FILES([tests/test125cellsPerturbed.sh], [chmod +x tests/test125cellsPerturbed.sh])
AC_CONFIG_FILES([tests/testPeriodicBC.sh], [chmod +x tests/testPeriodicBC.sh])
AC_CONFIG_FILES([tests/testPeriodicBCPerturbed.sh], [chmod +x tests/testPeriodicBCPerturbed.sh])
AC_CONFIG_FILES([tests/testInteractions.sh], [chmod +x tests/testInteractions.sh])
AC_CONFIG_FILES([tests/testParser.sh], [chmod +x tests/testParser.sh])
AC_CONFIG_FILES([tests/testSelectOutput.sh], [chmod +x tests/testSelectOutput.sh])
AC_CONFIG_FILES([tests/testFormat.sh], [chmod +x tests/testFormat.sh])

# Save the compilation options
AC_DEFINE_UNQUOTED([SWIFT_CONFIG_FLAGS],["$swift_config_flags"],[Flags passed to configure])

# Make sure the latest git revision string gets included, when we are
# working in a checked out repository.
test -d ${srcdir}/.git && touch ${srcdir}/src/version.c

#  Need to define this, instead of using fifth argument of AC_INIT, until
#  2.64. Defer until now as this redefines PACKAGE_URL, which can emit a
#  compilation error when testing with -Werror.
AC_DEFINE([PACKAGE_URL],["www.swiftsim.com"], [Package web pages])

# Generate output.
AC_OUTPUT

# Report general configuration.
AC_MSG_RESULT([
 ------- Summary --------

   $PACKAGE_NAME v.$PACKAGE_VERSION

   Compiler             : $CC
    - vendor            : $ax_cv_c_compiler_vendor
    - version           : $ax_cv_c_compiler_version
    - flags             : $CFLAGS
   MPI enabled          : $enable_mpi
   HDF5 enabled         : $with_hdf5
    - parallel          : $have_parallel_hdf5
   METIS/ParMETIS       : $have_metis / $have_parmetis
   FFTW3 enabled        : $have_fftw
   GSL enabled          : $have_gsl
   libNUMA enabled      : $have_numa
   GRACKLE enabled      : $have_grackle
   Special allocators   : $have_special_allocator
   CPU profiler         : $have_profiler
   Pthread barriers     : $have_pthread_barrier
   VELOCIraptor enabled : $have_velociraptor
   Particle Logger      : $with_logger

   Hydro scheme       : $with_hydro
   Dimensionality     : $with_dimension
   Kernel function    : $with_kernel
   Equation of state  : $with_eos
   Adiabatic index    : $with_gamma
   Riemann solver     : $with_riemann

   Gravity scheme      : $with_gravity
   Multipole order     : $with_multipole_order
   No gravity below ID : $no_gravity_below_id
   Make gravity glass  : $gravity_glass_making
   External potential  : $with_potential

   Entropy floor        : $with_entropy_floor
   Cooling function     : $with_cooling
   Chemistry            : $with_chemistry
   Tracers              : $with_tracers
   Stellar model        : $with_stars
   Star formation model : $with_star_formation
   Feedback model       : $with_feedback

   Individual timers           : $enable_timers
   Task debugging              : $enable_task_debugging
   Threadpool debugging        : $enable_threadpool_debugging
   Debugging checks            : $enable_debugging_checks
   Interaction debugging       : $enable_debug_interactions
   Stars interaction debugging : $enable_debug_interactions_stars
   Naive interactions          : $enable_naive_interactions
   Gravity checks              : $gravity_force_checks
   Custom icbrtf               : $enable_custom_icbrtf

 ------------------------])
