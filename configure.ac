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
AC_INIT([SWIFT],[1.0.0],[https://gitlab.cosma.dur.ac.uk/swift/swiftsim])
swift_config_flags="$*"

AC_COPYRIGHT
AC_CONFIG_SRCDIR([src/space.c])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([subdir-objects])

# Add local macro collection.
AC_CONFIG_MACRO_DIR([m4])

# Stop default CFLAGS from anyone except the environment.
: ${CFLAGS=""}

# Generate header file.
AC_CONFIG_HEADERS([config.h])


# Find and test the compiler.
AX_CHECK_ENABLE_DEBUG
# Enable POSIX and platform extension preprocessor macros.
AC_USE_SYSTEM_EXTENSIONS
AC_PROG_CC
AM_PROG_CC_C_O

# We need this for compilation hints and possibly FFTW.
AX_OPENMP

# If debug is selected then we also define SWIFT_DEVELOP_MODE to control
# any developer code options.
if test "x$ax_enable_debug" != "xno"; then
   AC_DEFINE([SWIFT_DEVELOP_MODE],1,[Enable developer code options])
fi

# C++ in GCC 6 and above has an issue with undefined the min() and max()
# macros. This hack works around that.
AC_DEFINE([_GLIBCXX_INCLUDE_NEXT_C_HEADERS],1,[Hack for min() and max() using g++ 6+])

# Enable POSIX and platform extension preprocessor macros.
AC_USE_SYSTEM_EXTENSIONS

# Check for C compiler version and vendor.
AX_COMPILER_VENDOR
AX_COMPILER_VERSION

# Check whether we have a recent enough GCC
if test "$ax_cv_c_compiler_vendor" = "gnu"; then
   AX_COMPARE_VERSION([$ax_cv_c_compiler_version], [ge], [8.1.0],
		      [gcc_handles_avx512="yes"],
		      [gcc_handles_avx512="no"])
fi

#  Restrict support.
AC_C_RESTRICT

# CSDS
AC_ARG_ENABLE([csds],
	[AS_HELP_STRING([--enable-csds],
		[enable the Continuous Simulation Data Stream]
	)],
	[with_csds="${enableval}"],
	[with_csds="no"]
)

if test "$with_csds" = "yes"; then
   AC_DEFINE([WITH_CSDS], 1, [csds enabled])
   # Ensure that the submodule is initialized
   ${srcdir}/tools/update-modules csds
   # The csds requires that long long is a 64bit type, let's
   # check that.
   AC_CHECK_SIZEOF([long long int])
   if test "$ac_cv_sizeof_long_long_int" != "8"; then
      AC_MSG_ERROR([The CSDS requires that 'long long int' has size 8 bytes])
   fi

   AC_CONFIG_SUBDIRS([csds])

fi
AM_CONDITIONAL([HAVECSDS],[test $with_csds = "yes"])


# Use best known optimization for the current architecture. Actual optimization
# happens later so we can avoid any issues it introduces with the compiler,
# but we need to know if that will happen now.
AC_ARG_ENABLE([optimization],
   [AS_HELP_STRING([--enable-optimization],
     [Enable compile time optimization flags for host @<:@yes/no@:>@]
   )],
   [enable_opt="$enableval"],
   [enable_opt="yes"]
)

# Interprocedural optimization support. Can need special handling for linking
# and archiving as well as compilation. Needs to be done before libtool is
# configured so we use the correct LD. It can give good improvements for
# clang based compilers, so the default is enabled, but we make that
# disabled when not optimizing or debugging support is enabled.
enable_ipo_default="yes";
if test "x$ax_enable_debug" = "xyes" -o "x$enable_opt" != "xyes"; then
   enable_ipo_default="no";
   AC_MSG_WARN([Interprocedural optimization support default is changed to false])
fi
AC_ARG_ENABLE([ipo],
   [AS_HELP_STRING([--enable-ipo],
     [Enable interprocedural optimization [default=yes unless debugging]]
   )],
   [enable_ipo="$enableval"],
   [enable_ipo="$enable_ipo_default"]
)

if test "$enable_ipo" = "yes"; then
   if test "$ax_cv_c_compiler_vendor" = "intel"; then
      CFLAGS="$CFLAGS -ip -ipo"
      LDFLAGS="$LDFLAGS -ipo"
      AC_CHECK_PROGS([AR], [xiar])
      AC_CHECK_PROGS([LD], [xild])
      AC_MSG_RESULT([added Intel interprocedural optimization support])
   elif test "$ax_cv_c_compiler_vendor" = "oneapi"; then
      CFLAGS="$CFLAGS -ipo"
      LDFLAGS="$LDFLAGS -ipo"
      AC_CHECK_PROGS([AR], [xiar])
      AC_CHECK_PROGS([LD], [xild])
      AC_CHECK_PROGS([RANLIB], [llvm-ranlib])
      AC_MSG_RESULT([added oneapi interprocedural optimization support])
   elif test "$ax_cv_c_compiler_vendor" = "gnu"; then
      CFLAGS="$CFLAGS -flto"
      LDFLAGS="$LDFLAGS -flto"
      AX_COMPARE_VERSION($ax_cv_c_compiler_version, [ge], [5.0.0],
         [
         AC_CHECK_PROGS([AR], [gcc-ar])
         AC_CHECK_PROGS([RANLIB], [gcc-ranlib])
         ], [:] )
      AC_MSG_RESULT([added GCC interprocedural optimization support])
   elif test "$ax_cv_c_compiler_vendor" = "clang"; then
      CFLAGS="$CFLAGS -flto"
      LDFLAGS="$LDFLAGS -flto"
      AC_CHECK_PROGS([RANLIB], [llvm-ranlib])
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
# Use extra flags set by AC_PROG_CC as part of $CC. Currently undocumented as ac_cv_prog_cc_stdc,
# so could change.
good_mpi="yes"
if test "$enable_mpi" = "yes"; then
    AX_MPI([CC="$MPICC $ac_cv_prog_cc_stdc" AC_DEFINE(HAVE_MPI, 1, [Define if you have the MPI library.]) ], [enable_mpi="no"])
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

# Add libtool support (now that CC is defined). Disable shared libraries by default.
LT_INIT([disable-shared])

# Need C99 and inline support. Only for autoconfs to version 2.69.
m4_version_prereq([2.71], [], [AC_PROG_CC_C99])

# Need inline support.
AC_C_INLINE


# If debugging try to show inlined functions.
if test "x$ax_enable_debug" = "xyes"; then
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
   elif test "$ax_cv_c_compiler_vendor" = "oneapi"; then
      CFLAGS="$CFLAGS -debug inline-debug-info"
   fi
fi

# Are we using the regular tasks (with atomics and no task conflicts) or
# are we using the atomic-free task-conflicting (slower) version?
AC_ARG_ENABLE([atomics-within-tasks],
   [AS_HELP_STRING([--disable-atomics-within-tasks],
     [Disable the use of atomic operations within tasks. This creates more task conflicts
      but allows for atomic-free and lock-free code within the tasks. This likely slows down
      the code @<:@yes/no@:>@]
   )],
   [enable_atomics_within_tasks="$enableval"],
   [enable_atomics_within_tasks="yes"]
)
AS_IF([test "x$enable_atomics_within_tasks" != "xno"],
   , # Note: atomics are allowed by default, we define the macro if we don't want them.
   [AC_DEFINE([SWIFT_TASKS_WITHOUT_ATOMICS],1,[Makes SWIFT use atomic-free and lock-free tasks.])
])


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
   LDFLAGS="$LDFLAGS -rdynamic -ldl"
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

# Check if cell graph is on.
AC_ARG_ENABLE([cell-graph],
   [AS_HELP_STRING([--enable-cell-graph],
     [Activate the cell graph @<:@yes/no@:>@]
   )],
   [enable_cell_graph="$enableval"],
   [enable_cell_graph="no"]
)
if test "$enable_cell_graph" = "yes"; then
   AC_DEFINE([SWIFT_CELL_GRAPH],1,[Enable cell graph])
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

# Check whether we want to default to naive cell interactions (stars)
AC_ARG_ENABLE([naive-interactions-stars],
   [AS_HELP_STRING([--enable-naive-interactions-stars],
     [Activate use of naive cell interaction functions for stars @<:@yes/no@:>@]
   )],
   [enable_naive_interactions_stars="$enableval"],
   [enable_naive_interactions_stars="no"]
)
if test "$enable_naive_interactions_stars" = "yes"; then
   AC_DEFINE([SWIFT_USE_NAIVE_INTERACTIONS_STARS],1,[Enable use of naive cell interaction functions for stars])
fi

# Check whether we want to default to naive cell interactions (stars)
AC_ARG_ENABLE([naive-interactions-rt],
   [AS_HELP_STRING([--enable-naive-interactions-rt],
     [Activate use of naive cell interaction functions for stars in RT tasks@<:@yes/no@:>@]
   )],
   [enable_naive_interactions_rt="$enableval"],
   [enable_naive_interactions_rt="no"]
)
if test "$enable_naive_interactions_rt" = "yes"; then
   AC_DEFINE([SWIFT_USE_NAIVE_INTERACTIONS_RT],1,[Enable use of naive cell interaction functions for stars in RT tasks])
fi

# Check if gravity force checks are on for some particles.
AC_ARG_ENABLE([gravity-force-checks],
   [AS_HELP_STRING([--enable-gravity-force-checks=<N>],
     [Activate expensive brute-force gravity checks for a fraction 1/N of all particles @<:@N@:>@]
   )],
   [gravity_force_checks="$enableval"],
   [gravity_force_checks="no"]
)
if test "$gravity_force_checks" = "yes"; then
   AC_MSG_ERROR(Need to specify the fraction of particles to check when using --enable-gravity-force-checks!)
elif test "$gravity_force_checks" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_GRAVITY_FORCE_CHECKS], [$enableval] ,[Enable gravity brute-force checks])
fi

# Check if hydro density checks are on for some particles.
AC_ARG_ENABLE([hydro-density-checks],
   [AS_HELP_STRING([--enable-hydro-density-checks],
     [Activate expensive brute-force hydro density checks for a fraction 1/N of all particles @<:@N@:>@]
   )],
   [hydro_density_checks="$enableval"],
   [hydro_density_checks="no"]
)
if test "$hydro_density_checks" = "yes"; then
   AC_MSG_ERROR(Need to specify the fraction of particles to check when using --enable-hydro-density-checks!)
elif test "$hydro_density_checks" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_HYDRO_DENSITY_CHECKS], [$enableval] ,[Enable hydro density brute-force checks])
fi

# Check if stars density checks are on for some particles.
AC_ARG_ENABLE([stars-density-checks],
   [AS_HELP_STRING([--enable-stars-density-checks],
     [Activate expensive brute-force stars density checks for a fraction 1/N of all particles @<:@N@:>@]
   )],
   [stars_density_checks="$enableval"],
   [stars_density_checks="no"]
)
if test "$stars_density_checks" = "yes"; then
   AC_MSG_ERROR(Need to specify the fraction of particles to check when using --enable-stars-density-checks!)
elif test "$stars_density_checks" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_STARS_DENSITY_CHECKS], [$enableval] ,[Enable stars density brute-force checks])
fi

# Check if ghost statistics are enabled
AC_ARG_ENABLE([ghost-statistics],
   [AS_HELP_STRING([--enable-ghost-statistics],
     [Gather statistics about the ghost iterations for hydro, stars and black holes in N bins @<:@N@:>@]
   )],
   [ghost_stats="$enableval"],
   [ghost_stats="no"]
)
if test "$ghost_stats" == "yes"; then
   AC_MSG_ERROR(Need to specify the number of bins when using --enable-ghost-statistics!)
elif test "$ghost_stats" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_GHOST_STATS], [$enableval] ,[Enable ghost statistics for hydro, stars and black holes])
fi

# Check whether we want to switch on glass making
AC_ARG_ENABLE([glass-making],
   [AS_HELP_STRING([--enable-glass-making],
     [Activate the glass-making procedure by reversing the sign of gravity @<:@yes/no@:>@]
   )],
   [gravity_glass_making="$enableval"],
   [gravity_glass_making="no"]
)
if test "$gravity_glass_making" = "yes"; then
   AC_DEFINE([SWIFT_MAKE_GRAVITY_GLASS], 1, [Make the code run in a way to produce a glass file for gravity/cosmology])
fi

# Check if we want to zero the gravity forces for all particles below some ID.
AC_ARG_ENABLE([no-gravity-below-id],
   [AS_HELP_STRING([--enable-no-gravity-below-id=<N>],
     [Zeros the gravitational acceleration of all particles with an ID smaller than @<:@N@:>@]
   )],
   [no_gravity_below_id="$enableval"],
   [no_gravity_below_id="no"]
)
if test "$no_gravity_below_id" = "yes"; then
   AC_MSG_ERROR(Need to specify the ID below which particles get zero forces when using --enable-no-gravity-below-id!)
elif test "$no_gravity_below_id" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_NO_GRAVITY_BELOW_ID], [$enableval] ,[Particles with smaller ID than this will have zero gravity forces])
fi

# Check if we want to use boundary particles.
AC_ARG_ENABLE([boundary-particles],
   [AS_HELP_STRING([--enable-boundary-particles=<N>],
     [Set all particles with an ID smaller than @<:@N@:>@ as boundary particles (i.e. receive zero gravity + hydro forces).]
   )],
   [boundary_particles="$enableval"],
   [boundary_particles="no"]
)
if test "$boundary_particles" = "yes"; then
   AC_MSG_ERROR(Need to specify the ID below which particles get zero forces when using --enable-boundary-particles!)
elif test "$boundary_particles" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_NO_GRAVITY_BELOW_ID], [$enableval] ,[Particles with smaller ID than this will have zero gravity forces])
   AC_DEFINE_UNQUOTED([SWIFT_BOUNDARY_PARTICLES], [$enableval] ,[Particles with smaller ID than this will be considered as boundaries.])
fi

# Check if we want to use fixed boundary particles.
AC_ARG_ENABLE([fixed-boundary-particles],
   [AS_HELP_STRING([--enable-fixed-boundary-particles=<N>],
     [Set all particles with an ID smaller than @<:@N@:>@ as fixed boundary particles (i.e. receive zero gravity + hydro forces + zero velocity), this mode enables also --enable-boundary-particles and --enable-no-gravity-below-id.]
   )],
   [fixed_boundary_particles="$enableval"],
   [fixed_boundary_particles="no"]
)
if test "$fixed_boundary_particles" = "yes"; then
   AC_MSG_ERROR(Need to specify the ID below which particles get zero forces when using --enable-boundary-particles!)
elif test "$fixed_boundary_particles" != "no"; then
   AC_DEFINE_UNQUOTED([SWIFT_NO_GRAVITY_BELOW_ID], [$enableval] ,[Particles with smaller ID than this will have zero gravity forces])
   AC_DEFINE_UNQUOTED([SWIFT_BOUNDARY_PARTICLES], [$enableval] ,[Particles with smaller ID than this will be considered as boundaries.])
   AC_DEFINE_UNQUOTED([SWIFT_FIXED_BOUNDARY_PARTICLES], [$enableval] ,[Particles with smaller ID than this will be considered as boundaries.])
fi

# Check if fixed entropy is on for settling planetary initial conditions
AC_ARG_ENABLE([planetary-fixed-entropy],
   [AS_HELP_STRING([--enable-planetary-fixed-entropy],
     [Force entropies to stay fixed for settling planetary initial conditions @<:@yes/no@:>@]
   )],
   [planetary_fixed_entropy="$enableval"],
   [planetary_fixed_entropy="no"]
)
if test "$planetary_fixed_entropy" = "yes"; then
   AC_DEFINE([PLANETARY_FIXED_ENTROPY],1,[Enable planetary fixed entropy])
fi

# Check whether we have any of the ARM v8.1 tick timers
AX_ASM_ARM_PMCCNTR
AX_ASM_ARM_CNTVCT

# See if we want memuse reporting.
AC_ARG_ENABLE([memuse-reports],
   [AS_HELP_STRING([--enable-memuse-reports],
     [Output reports about significant memory allocations@<:@yes/no@:>@]
   )],
   [enable_memuse_reports="$enableval"],
   [enable_memuse_reports="no"]
)
if test "$enable_memuse_reports" = "yes"; then
   AC_DEFINE([SWIFT_MEMUSE_REPORTS],1,[Enable memory usage reports])
fi

# The system memory report depends on the presence of /proc/self/statm, i.e.
# a linux OS, check for that.
if test -f "/proc/self/statm"; then
   AC_DEFINE([SWIFT_MEMUSE_STATM],1,[Have /proc/self/statm capability])
fi

# Check if we want to make the dumper thread active.
AC_ARG_ENABLE([dumper],
   [AS_HELP_STRING([--enable-dumper],
     [Dump active tasks and memory use (if configured)@<:@yes/no@:>@]
   )],
   [enable_dumper="$enableval"],
   [enable_dumper="no"]
)
if test "$enable_dumper" = "yes"; then
   AC_DEFINE([SWIFT_DUMPER_THREAD],1,[Enable dumper thread])
fi

# See if we want mpi reporting.
AC_ARG_ENABLE([mpiuse-reports],
   [AS_HELP_STRING([--enable-mpiuse-reports],
     [Output reports about MPI tasks requests@<:@yes/no@:>@]
   )],
   [enable_mpiuse_reports="$enableval"],
   [enable_mpiuse_reports="no"]
)
if test "$enable_mpiuse_reports" = "yes"; then
   AC_DEFINE([SWIFT_MPIUSE_REPORTS],1,[Enable MPI task reports])
fi


# Define HAVE_POSIX_MEMALIGN if it works.
AX_FUNC_POSIX_MEMALIGN

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

# Only optimize if allowed, otherwise assume user will set CFLAGS as
# appropriate. Note argument check is done earlier so we can configure
# other options related to optimization.
if test "$enable_opt" = "yes" ; then

   # Choose the best flags for this compiler and architecture
   ac_test_CFLAGS="no"
   AX_CC_MAXOPT
   ac_test_CFLAGS="yes"

   # Choose the best flags for the gravity sub-library on this compiler and
   # architecture. Note we use OpenMP as a compiler hints for loop vectorization.
   GRAVITY_CFLAGS="$GRAVITY_CFLAGS $OPENMP_CFLAGS"
   if test "$ax_cv_c_compiler_vendor" = "intel"; then
      case "$icc_flags" in
      	 *CORE-AVX512*)
            GRAVITY_CFLAGS="$GRAVITY_CFLAGS -qopt-zmm-usage=high"
	    ;;
	 *)
	    AC_MSG_NOTICE([No additional flags needed for gravity on this platform])
	    ;;
      esac
   elif test "$ax_cv_c_compiler_vendor" = "oneapi"; then
      case "$icc_flags" in
         *CORE-AVX512*)
            GRAVITY_CFLAGS="$GRAVITY_CFLAGS -qopt-zmm-usage=high"
            ;;
         *)
            AC_MSG_NOTICE([No additional flags needed for gravity on this platform])
            ;;
      esac
   elif test "$ax_cv_c_compiler_vendor" = "gnu"; then
      if test "$gcc_handles_avx512" = "yes"; then
         case "$ax_cv_gcc_archflag" in
	    *skylake-avx512*)
               GRAVITY_CFLAGS="$GRAVITY_CFLAGS -mprefer-vector-width=512"
	       ;;
	    *)
	       AC_MSG_NOTICE([No additional flags needed for gravity on this platform])
	       ;;
         esac
      else
         AC_MSG_NOTICE([No additional flags needed for gravity on this platform])
      fi
   elif test "$ax_cv_c_compiler_vendor" = "clang"; then
      #  Could be a number of compilers. Check for aocc specific flags we want
      #  to use.
      AX_CHECK_COMPILE_FLAG("-zopt", [GRAVITY_CFLAGS="$GRAVITY_CFLAGS -fvectorize -zopt"])
      case "$ax_cv_gcc_archflag" in
	    *skylake-avx512*)
               GRAVITY_CFLAGS="$GRAVITY_CFLAGS -mprefer-vector-width=512"
	       ;;
	    *znver[[4-9]])
               GRAVITY_CFLAGS="$GRAVITY_CFLAGS -mprefer-vector-width=512"
               ;;
            *)
               :
            ;;
       esac
   else
      AC_MSG_WARN([Do not know what best gravity vectorization flags to choose for this compiler])
   fi
   AC_ARG_VAR([GRAVITY_CFLAGS], [C compiler flags added to the basic CFLAGS to compile
   				 the gravity-related files.])

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
      elif test "$ax_cv_c_compiler_vendor" = "oneapi"; then
         CFLAGS="$CFLAGS -no-vec"
         AC_MSG_RESULT([disabled oneAPI vectorization])
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
      AC_MSG_RESULT([Adding address sanitizer support... yes])

        # Check if we have access to the __lsan_ignore_object() call for
        # marking memory allocations as deliberately leaked.
        AC_LINK_IFELSE([AC_LANG_SOURCE([[
                           #include <stdlib.h>
                           #include <sanitizer/lsan_interface.h>
                           int main(int argc, char *argv[]) {
                              void *p = malloc(1);
                              __lsan_ignore_object(p);
                              return 0;
                           }]])],
           [AC_DEFINE(HAVE_LSAN_IGNORE_OBJECT, 1, [Have __lsan_ignore_object() call])],
           [AC_MSG_WARN([Sanitizer enabled but no __lsan_ignore_object()])])
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


# MPI mesh gravity
AC_ARG_ENABLE([mpi-mesh-gravity],
	[AS_HELP_STRING([--enable-mpi-mesh-gravity],
		[enable parallel mesh gravity (requires FFTW MPI library) @<:@no/yes@:>@]
	)],
	[with_mpi_mesh_gravity="${enableval}"],
	[with_mpi_mesh_gravity="no"]
)
# Autoconf stuff.
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_EGREP


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
if test "x$have_pthread_barrier" = "xyes"; then
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

   #  Recent METIS releases have an external GKlib, test for that before
   # giving up. Assume sane and in the same directories.
   if test "x$have_metis" = "xno"; then
      if test "x$with_metis" != "xyes" -a "x$with_metis" != "x"; then
         METIS_LIBS="-L$with_metis/lib -lmetis -lGKlib"
         METIS_INCS="-I$with_metis/include"
      else
         METIS_LIBS="-lmetis -lGKlib"
         METIS_INCS=""
      fi
      AC_CHECK_LIB([metis],[METIS_PartGraphRecursive], [have_metis="yes"],
                   [have_metis="no"], $METIS_LIBS)
   fi

   if test "$have_metis" = "yes"; then
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

# A build may use an external METIS library, check for that.
   if test "$have_parmetis" = "no"; then
      if test "x$with_parmetis" != "xyes" -a "x$with_parmetis" != "x"; then
         PARMETIS_LIBS="-L$with_parmetis/lib -lparmetis -lmetis"
         PARMETIS_INCS="-I$with_parmetis/include"
      else
         PARMETIS_LIBS="-lparmetis -lmetis"
         PARMETIS_INCS=""
      fi
      # Note use different function to avoid caching of first check.
      AC_CHECK_LIB([parmetis],[ParMETIS_V3_PartKway], [have_parmetis="yes"],
                   [have_parmetis="no"], [$METIS_LIBS $PARMETIS_LIBS])
   fi

# A build may use an external GKlib in later releases...
   if test "$have_parmetis" = "no"; then
      if test "x$with_parmetis" != "xyes" -a "x$with_parmetis" != "x"; then
         PARMETIS_LIBS="-L$with_parmetis/lib -lparmetis -lGKlib"
         PARMETIS_INCS="-I$with_parmetis/include"
      else
         PARMETIS_LIBS="-lparmetis -lGKlib"
         PARMETIS_INCS=""
      fi
      # Note use different function to avoid caching of first check.
      AC_CHECK_LIB([parmetis],[ParMETIS_V3_PartGeom], [have_parmetis="yes"],
                   [have_parmetis="no"], [$METIS_LIBS $PARMETIS_LIBS])

   fi
   if test "$have_parmetis" = "yes"; then
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
# If FFTW is found, we check whether this is the threaded or openmp version.
have_fftw="no"
have_mpi_fftw="no"
have_threaded_fftw="no"
have_openmp_fftw="no"
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

      else

         # Same checks for OpenMP if preferred threaded failed.
         if test "x$with_fftw" != "xyes" -a "x$with_fftw" != "xtest" -a "x$with_fftw" != "x"; then
            FFTW_OPENMP_LIBS="-L$with_fftw/lib -lfftw3_omp -lfftw3"
            FFTW_OPENMP_INCS="-I$with_fftw/include"
         else
            FFTW_OPENMP_LIBS="-lfftw3_omp -lfftw3"
            FFTW_OPENMP_INCS=""
         fi

         # Verify that the library works. Note requires AX_OPENMP called above.
         old_CFLAGS=$CFLAGS
         CFLAGS="$CFLAGS $OPENMP_CFLAGS"
         AC_CHECK_LIB([fftw3],[fftw_init_threads],[have_openmp_fftw="yes"],
		      [have_openmp_fftw="no"], $FFTW_OPENMP_LIBS)

         # If found, update things
         if test "x$have_openmp_fftw" = "xyes"; then
            # Note OpenMP and pthreads use mostly the same calls, so define both.
            AC_DEFINE([HAVE_THREADED_FFTW],1,[The threaded OpenMP FFTW library appears to be present.])
            AC_DEFINE([HAVE_OPENMP_FFTW],1,[The OpenMP FFTW library appears to be present.])
            FFTW_LIBS=$FFTW_OPENMP_LIBS
            FFTW_INCS=$FFTW_OPENMP_INCS
         else
            # Put CFLAGS back.
            CFLAGS=$old_CFLAGS
         fi
      fi

   fi

   # If MPI mesh gravity is not disabled, check whether we have the MPI version of FFTW
   if test "x$enable_mpi" = "xyes" -a "x$with_mpi_mesh_gravity" != "xno"; then
      # Was FFTW's location specifically given?
      if test "x$with_fftw" != "xyes" -a "x$with_fftw" != "xtest" -a "x$with_fftw" != "x"; then
         FFTW_MPI_LIBS="-L$with_fftw/lib -lfftw3_mpi -lfftw3"
         FFTW_MPI_INCS="-I$with_fftw/include"
      else
         FFTW_MPI_LIBS="-lfftw3_mpi -lfftw3"
         FFTW_MPI_INCS=""
      fi

      # Verify that the library has MPI support
      AC_CHECK_LIB([fftw3],[fftw_mpi_init],[have_mpi_fftw="yes"],
                   [have_mpi_fftw="no"], $FFTW_MPI_LIBS)

      # If found, update things. Don't add MPI flags to FFTW_*_LIBS etc because
      # we don't want to link the MPI library into the non-MPI swift executable.
      if test "x$have_mpi_fftw" = "xyes"; then
         AC_DEFINE([HAVE_MPI_FFTW],1,[The MPI FFTW library appears to be present.])
      else
         if test "x$with_mpi_mesh_gravity" = "xyes" ; then
            AC_MSG_ERROR("Unable to find FFTW MPI library for MPI mesh gravity")
         fi
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
have_arm_fftw="no"
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
      fi
   # FFTW was specified, check that it was a valid location.
   else
      AC_CHECK_LIB([armpl_lp64],[fftw_malloc],
         AC_DEFINE([HAVE_FFTW],1,[The FFTW library appears to be present.]),
         AC_MSG_ERROR(something is wrong with the FFTW library!), $FFTW_LIBS)
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

AC_SUBST([FFTW_MPI_LIBS])
AC_SUBST([FFTW_MPI_INCS])
AM_CONDITIONAL([HAVEMPIFFTW],[test -n "$FFTW_MPI_LIBS"])

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
        intel | gnu | clang | oneapi)
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
        intel | gnu | clang | oneapi)
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
        intel | gnu | clang | oneapi)
             CFLAGS="$CFLAGS -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free"
          ;;
      esac

   else
      TBBMALLOC_LIBS=""
   fi
fi
AC_SUBST([TBBMALLOC_LIBS])
AM_CONDITIONAL([HAVETBBMALLOC],[test -n "$TBBMALLOC_LIBS"])

# check for a random seed
AC_ARG_WITH([random-seed],
    [AS_HELP_STRING([--with-random-seed=SHORT INT],
       [Set the random seed.]
    )],
    [with_random_seed="$withval"],
    [with_random_seed="0"]
)
AC_DEFINE_UNQUOTED([SWIFT_RANDOM_SEED_XOR], [$with_random_seed],[Value of the random seed.])


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

   AS_VAR_APPEND([GRACKLE_LIBS], ["$FCLIBS"])

   AC_CHECK_LIB(
      [grackle],
      [initialize_chemistry_data],
      [AC_DEFINE([HAVE_GRACKLE],1,[The GRACKLE library appears to be present.])
        AC_DEFINE([CONFIG_BFLOAT_8],1,[Use doubles in grackle])
      ],
      [AC_MSG_ERROR(Cannot find grackle library! Please consult the documentation for specific version required.)],
      [$GRACKLE_LIBS])

   AC_CHECK_LIB(
      [grackle],
      [set_velocity_units],
      [ : ], # : means do nothing. Leaving this argument empty triggers AC default behaviour, which breaks stuff down the line.
      [AC_MSG_ERROR(Wrong grackle library version. Please consult the documentation for specifics.)],
      [$GRACKLE_LIBS])

   AC_CHECK_LIB(
      [grackle],
      [get_grackle_version],
      [ : ], # : means do nothing
      [AC_MSG_ERROR(Wrong grackle library version. Please consult the documentation for specifics.)],
      [$GRACKLE_LIBS])


fi
AC_SUBST([GRACKLE_LIBS])
AC_SUBST([GRACKLE_INCS])
AM_CONDITIONAL([HAVEGRACKLE],[test -n "$GRACKLE_LIBS"])

# Check for VELOCIraptor, non-MPI.
have_velociraptor="no"
AC_ARG_WITH([velociraptor],
    [AS_HELP_STRING([--with-velociraptor=PATH],
       [Directory where velociraptor library exists @<:@yes/no@:>@]
    )],
    [with_velociraptor="$withval"],
    [with_velociraptor="no"]
)
if test "x$with_velociraptor" != "xno"; then
   if test "x$with_velociraptor" != "xyes" -a "x$with_velociraptor" != "x"; then
      VELOCIRAPTOR_LIBS="-L$with_velociraptor -lvelociraptor -lstdc++ -lhdf5"
      CFLAGS="$CFLAGS -fopenmp"
   else
      VELOCIRAPTOR_LIBS=""
   fi

   have_velociraptor="yes"
   AS_VAR_APPEND([VELOCIRAPTOR_LIBS], ["$FCLIBS"])

   AC_CHECK_LIB(
      [velociraptor],
      [InitVelociraptor],
      [AC_DEFINE([HAVE_VELOCIRAPTOR],1,[The non-MPI VELOCIraptor library appears to be present.])],
      [AC_MSG_ERROR(Cannot find non-MPI VELOCIraptor library at $with_velociraptor or incompatible HDF5 library loaded.)],
      [$VELOCIRAPTOR_LIBS $HDF5_LDFLAGS $HDF5_LIBS $GSL_LIBS]
   )
fi
AC_SUBST([VELOCIRAPTOR_LIBS])
AM_CONDITIONAL([HAVEVELOCIRAPTOR],[test -n "$VELOCIRAPTOR_LIBS"])

# Now that we found VELOCIraptor, let's check how it was compiled.
if test "$have_velociraptor" = "yes"; then
    AC_CHECK_LIB(
       [velociraptor],
       [VR_NOMASS],
       [AC_DEFINE([HAVE_VELOCIRAPTOR_WITH_NOMASS],1,[The VELOCIraptor library has been compiled with the NOMASS option. Only useful if running a uniform box.])],
       [AC_MSG_RESULT(VELOCIraptor not compiled to so as to *not* store masses per particle.)],
       [$VELOCIRAPTOR_LIBS $HDF5_LDFLAGS $HDF5_LIBS $GSL_LIBS]
    )
fi

#  Check for MPI VELOCIraptor, same as above.
have_mpi_velociraptor="no"
AC_ARG_WITH([velociraptor-mpi],
    [AS_HELP_STRING([--with-velociraptor-mpi=PATH],
       [Directory where MPI version of velociraptor library exists @<:@yes/no@:>@]
    )],
    [with_mpi_velociraptor="$withval"],
    [with_mpi_velociraptor="no"]
)
if test "x$with_mpi_velociraptor" != "xno"; then
   if test "x$with_mpi_velociraptor" != "xyes" -a "x$with_mpi_velociraptor" != "x"; then
      VELOCIRAPTOR_MPI_LIBS="-L$with_mpi_velociraptor -lvelociraptor -lmpi -lstdc++ -lhdf5"
      CFLAGS="$CFLAGS -fopenmp"
   else
      VELOCIRAPTOR_MPI_LIBS=""
   fi

   have_mpi_velociraptor="yes"
   AS_VAR_APPEND([VELOCIRAPTOR_MPI_LIBS], ["$FCLIBS"])

   AC_CHECK_LIB(
      [velociraptor],
      [InitVelociraptor],
      [AC_DEFINE([HAVE_MPI_VELOCIRAPTOR],1,[The MPI VELOCIraptor library appears to be present.])],
      [AC_MSG_ERROR(Cannot find MPI VELOCIraptor library at $with_mpi_velociraptor or incompatible HDF5 library loaded.)],
      [$VELOCIRAPTOR_MPI_LIBS $HDF5_LDFLAGS $HDF5_LIBS $GSL_LIBS]
   )
fi
AC_SUBST([VELOCIRAPTOR_MPI_LIBS])
AM_CONDITIONAL([HAVEVELOCIRAPTOR],[test -n "$VELOCIRAPTOR_MPI_LIBS"])

# Let's check how this one was compiled.
if test "$have_mpi_velociraptor" = "yes"; then
    AC_CHECK_LIB(
       [velociraptor],
       [VR_NOMASS],
       [AC_DEFINE([HAVE_VELOCIRAPTOR_WITH_NOMASS],1,[The MPI VELOCIraptor library has been compiled with the NOMASS option. Only useful if running a uniform box.])],
       [AC_MSG_RESULT(VELOCIraptor not compiled to so as to *not* store masses per particle.)],
       [$VELOCIRAPTOR_MPI_LIBS $HDF5_LDFLAGS $HDF5_LIBS $GSL_LIBS]
    )
fi

# If we have one library, but not the other then use that for both.
if test "$have_mpi_velociraptor" = "yes" -a "$have_velociraptor" != "yes"; then
   VELOCIRAPTOR_LIBS="$VELOCIRAPTOR_MPI_LIBS"
   AC_SUBST([VELOCIRAPTOR_LIBS])
elif test "$have_velociraptor" = "yes" -a "$have_mpi_velociraptor" != "yes"; then
   VELOCIRAPTOR_MPI_LIBS="$VELOCIRAPTOR_LIBS"
   AC_SUBST([VELOCIRAPTOR_MPI_LIBS])
fi

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

# Check if we should be writing out most bound "orphan" particles from Velociraptor
AC_ARG_ENABLE([velociraptor-orphans],
    [AS_HELP_STRING([--enable-velociraptor-orphans],
       [Enable output of orphan particles @<:@yes/no@:>@]
    )],
    [enable_velociraptor_orphans="$enableval"],
    [enable_velociraptor_orphans="no"]
)
if test "$enable_velociraptor_orphans" = "yes"; then
   AC_DEFINE([HAVE_VELOCIRAPTOR_ORPHANS], 1, [Orphan particles should be written out])
fi

# Check if lightcone output is on.
AC_ARG_ENABLE([lightcone],
   [AS_HELP_STRING([--enable-lightcone],
     [Activate lightcone outputs.],
   )],
   [enable_lightcone="$enableval"],
   [enable_lightcone="no"]
)
if test "$enable_lightcone" = "yes"; then
   # Check for healpix for lightcone maps. May require cfitsio
   # This sets CHEALPIX_LIBS and CHEALPIX_CFLAGS and #defines HAVE_CHEALPIX.
   # It also adds a --with-cfitsio flag in case cfitsio is installed in a
   # different location from healpix.
   GV_FIND_LIBRARY([cfitsio], [CFITSIO], [cfitsio], [cfitsio], [ffclos])
   TMP_LIBS=${LIBS}
   LIBS="${CFITSIO_LIBS} ${LIBS}"
   GV_FIND_LIBRARY([chealpix], [CHEALPIX], [chealpix], [chealpix], [ang2vec])
   LIBS=${TMP_LIBS}
   have_chealpix=${USE_CHEALPIX}
   CHEALPIX_LIBS="${CHEALPIX_LIBS} ${CFITSIO_LIBS}"
   AC_DEFINE([WITH_LIGHTCONE], 1, [Enable lightcone outputs])
   if test "$have_chealpix" != "yes"; then
      AC_MSG_ERROR([Lightcone output requires the HEALPix C API. Please configure with --with-chealpix.])
   fi
   # Also need to make sure we have GSL if we're making lightcones
   if test "$have_gsl" != "yes"; then
      AC_MSG_ERROR([Lightcone output requires GSL. Please configure with --with-gsl.])
   fi
else
   have_chealpix="no"
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



# Check for Sundials (required for the SPHM1RT library).
# There is a problems with the headers of this library
# as they do not pass the strict prototypes check when
# installed outside of the system directories. So we
# need to do this check in two phases.
have_sundials="no"
SUNDIALS_LIBS=""
SUNDIALS_INCS=""
AC_ARG_WITH([sundials],
    [AS_HELP_STRING([--with-sundials=PATH],
       [root directory where sundials is installed @<:@yes/no@:>@]
    )],
    [with_sundials="$withval"],
    [with_sundials="no"]
)
if test "x$with_sundials" != "xno"; then
   AC_PROG_FC
   AC_FC_LIBRARY_LDFLAGS
   if test "x$with_sundials" != "xyes" -a "x$with_sundials" != "x"; then
      SUNDIALS_LIBS="-L$with_sundials/lib -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense"
      SUNDIALS_INCS="-I$with_sundials/include"
   else
      SUNDIALS_LIBS="-lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense"
      SUNDIALS_INCS=""
   fi

   AC_CHECK_LIB([sundials_cvode], [CVode], [have_sundials="yes"],
                [have_sundials="no"], $SUNDIALS_LIBS)

   if test "$have_sundials" == "yes"; then
      AC_DEFINE([HAVE_SUNDIALS],1,[The SUNDIALS library is present.])
   else
      if test "x$with_sundials" != "xyes" -a "x$with_sundials" != "x"; then
      	 # It might be that the libraries are in
      	 # /lib64 rather than /lib
      	 SUNDIALS_LIBS="-L$with_sundials/lib64 -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense"

	 # unset cached result of previous AC_CHECK_LIB
	 unset ac_cv_lib_sundials_cvode_CVode

   	 AC_CHECK_LIB([sundials_cvode], [CVode], [have_sundials="yes"], [have_sundials="no"], $SUNDIALS_LIBS)

   	 if test "$have_sundials" == "yes"; then
      	    AC_DEFINE([HAVE_SUNDIALS],1,[The SUNDIALS library is present.])
	 else
	    AC_MSG_ERROR("Failed to find a SUNDIALS library")
	 fi

      else
      	 AC_MSG_ERROR("Failed to find a SUNDIALS library")
      fi
   fi
fi
AC_SUBST([SUNDIALS_LIBS])


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
AC_CHECK_HEADERS_ONCE([sys/time.h])
if test $ac_cv_header_sys_time_h = yes; then
  AC_DEFINE([TIME_WITH_SYS_TIME],[1],[Define to 1 if you can safely include both <sys/time.h>
	     and <time.h>.  This macro is obsolete.])
fi

AC_CHECK_HEADERS([sys/time.h], [], [],
[#ifdef HAVE_SYS_TIME_H
# include <sys/time.h>
#endif
])
AC_CHECK_HEADERS([c_asm.h], [], [],
[#ifdef HAVE_C_ASM_H
# include <c_asm.h>
#endif
])
AC_CHECK_HEADERS([intrinsics.h], [], [],
[#ifdef HAVE_INTRINSICS_H
# include <intrinsics.h>
#endif
])
AC_CHECK_HEADERS([mach/mach_time.h], [], [],
[#ifdef HAVE_MACH_MACH_TIME_H
# include <mach/mach_time.h>
#endif
])

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
AC_ARG_ENABLE(armv7a-cntvct, [AS_HELP_STRING([--enable-armv7a-cntvct],[enable the cycle counter on Armv7a via the CNTVCT register])], have_armv7acntvct=$enableval)
if test "$have_armv7acntvct"x = "yes"x; then
      AC_DEFINE(HAVE_ARMV7A_CNTVCT,1,[Define if you have enabled the CNTVCT cycle counter on ARMv7a])
fi

AC_ARG_ENABLE(armv7a-pmccntr, [AS_HELP_STRING([--enable-armv7a-pmccntr],[enable the cycle counter on Armv7a via the PMCCNTR register])], have_armv7apmccntr=$enableval)
if test "$have_armv7apmccntr"x = "yes"x; then
      AC_DEFINE(HAVE_ARMV7A_PMCCNTR,1,[Define if you have enabled the PMCCNTR cycle counter on ARMv7a])
fi

# Check if we have native exp10 and exp10f functions. If not failback to our
# implementations. On Apple/CLANG we have __exp10, so also check for that
# if the compiler is clang.
AC_CHECK_LIB([m],[exp10], [AC_DEFINE([HAVE_EXP10],1,[The exp10 function is present.])])
AC_CHECK_LIB([m],[exp10f], [AC_DEFINE([HAVE_EXP10F],1,[The exp10f function is present.])])
if test "$ax_cv_c_compiler_vendor" = "clang"; then
      AC_CHECK_LIB([m],[__exp10], [AC_DEFINE([HAVE___EXP10],1,[The __exp10 function is present.])])
      AC_CHECK_LIB([m],[__exp10f], [AC_DEFINE([HAVE___EXP10F],1,[The __exp10f function is present.])])
fi

# Check if we have native sincos and sincosf functions. If not failback to our
# implementations. On Apple/CLANG we have __sincos, so also check for that
# if the compiler is clang.
AC_CHECK_LIB([m],[sincos], [AC_DEFINE([HAVE_SINCOS],1,[The sincos function is present.])])
AC_CHECK_LIB([m],[sincosf], [AC_DEFINE([HAVE_SINCOSF],1,[The sincosf function is present.])])
if test "$ax_cv_c_compiler_vendor" = "clang"; then
      AC_CHECK_LIB([m],[__sincos], [AC_DEFINE([HAVE___SINCOS],1,[The __sincos function is present.])])
      AC_CHECK_LIB([m],[__sincosf], [AC_DEFINE([HAVE___SINCOSF],1,[The __sincosf function is present.])])
fi
 
# The aocc compiler has optimized maths libraries that we should use. Check
# any clang for this support. Note do this after the basic check for maths
# as we need to make sure -lm follows. Also note needs -Ofast or -ffast-math
# so only when optimizing.
if test "$enable_opt" = "yes" -a "$ax_cv_c_compiler_vendor" = "clang"; then
   have_almfast="yes"
   AC_CHECK_LIB([almfast],[amd_fastexp],[LIBS="-fveclib=AMDLIBM -fsclrlib=AMDLIBM -lalmfast -lamdlibm $LIBS"],[have_almfast="no"],[-lamdlibm -lm])
   if test "$have_almfast" = "no"; then
      # Less optimized version.
      AC_CHECK_LIB([amdlibm],[sqrt],,,[-lm])
   fi
fi

# Check for glibc extension backtrace().
AC_CHECK_FUNCS([backtrace backtrace_symbols])

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
          gnu | clang | oneapi)
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
          intel | gnu | clang | oneapi)
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


# Second part of the Sundials library checks.
# We now decide if we need to use -isystem to
# get around the strict-prototypes problem. Assumes
# isystem is available when strict-prototypes is.
if test "x$with_sundials" != "xno"; then
   if test "x$with_sundials" != "xyes" -a "x$with_sundials" != "x"; then
        case "$CFLAGS" in
            *strict-prototypes*)
	        SUNDIALS_INCS="-isystem$with_sundials/include"
            ;;
            *)
                SUNDIALS_INCS="-I$with_sundials/include"
            ;;
        esac
   fi
fi
AC_SUBST([SUNDIALS_INCS])


# Various package configuration options.

# Master subgrid options
# If you add a restriction (e.g. no cooling, chemistry or hydro)
# you will need to check for overwrite after reading the additional options.
# As an example for this, see the call to AC_ARG_WITH for cooling.
AC_ARG_WITH([subgrid],
	[AS_HELP_STRING([--with-subgrid=<subgrid>],
		[Master switch for subgrid methods. Inexperienced user should start here. Options are: @<:@none, GEAR, AGORA, QLA, QLA-EAGLE, EAGLE, EAGLE-XL, SPIN_JET_EAGLE default: none@:>@]
	)],
	[with_subgrid="$withval"],
	[with_subgrid=none]
)

# Default values
with_subgrid_cooling=none
with_subgrid_chemistry=none
with_subgrid_tracers=none
with_subgrid_entropy_floor=none
with_subgrid_pressure_floor=none
with_subgrid_stars=none
with_subgrid_star_formation=none
with_subgrid_feedback=none
with_subgrid_sink=none
with_subgrid_extra_io=none

case "$with_subgrid" in
   yes)
      AC_MSG_ERROR([Invalid option. A subgrid model must be chosen.])
   ;;
   none)
   ;;
   GEAR)
	with_subgrid_cooling=grackle_0
	with_subgrid_chemistry=GEAR_10
	with_subgrid_pressure_floor=none
	with_subgrid_stars=GEAR
	with_subgrid_star_formation=GEAR
	with_subgrid_feedback=GEAR
	with_subgrid_black_holes=none
	with_subgrid_sink=none
	with_subgrid_extra_io=none
	enable_fof=no
   ;;
   AGORA)
	with_subgrid_cooling=grackle_0
	with_subgrid_chemistry=AGORA
	with_subgrid_pressure_floor=GEAR
	with_subgrid_stars=GEAR
	with_subgrid_star_formation=GEAR
	with_subgrid_feedback=AGORA
	with_subgrid_black_holes=none
	with_subgrid_sink=none
	with_subgrid_extra_io=none
	enable_fof=no
   ;;
   QLA)
	with_subgrid_cooling=QLA
	with_subgrid_chemistry=QLA
	with_subgrid_tracers=none
	with_subgrid_entropy_floor=QLA
	with_subgrid_stars=basic
	with_subgrid_star_formation=QLA
	with_subgrid_feedback=none
	with_subgrid_black_holes=none
	with_subgrid_sink=none
	with_subgrid_extra_io=none
	enable_fof=no
   ;;
   QLA-EAGLE)
	with_subgrid_cooling=QLA-EAGLE
	with_subgrid_chemistry=QLA
	with_subgrid_tracers=none
	with_subgrid_entropy_floor=QLA
	with_subgrid_stars=basic
	with_subgrid_star_formation=QLA
	with_subgrid_feedback=none
	with_subgrid_black_holes=none
	with_subgrid_sink=none
	enable_fof=no
   ;;
   EAGLE)
	with_subgrid_cooling=EAGLE
	with_subgrid_chemistry=EAGLE
	with_subgrid_tracers=EAGLE
	with_subgrid_entropy_floor=EAGLE
	with_subgrid_stars=EAGLE
	with_subgrid_star_formation=EAGLE
	with_subgrid_feedback=EAGLE
	with_subgrid_black_holes=EAGLE
	with_subgrid_sink=none
	with_subgrid_extra_io=none
	enable_fof=yes
   ;;
   EAGLE-XL)
	with_subgrid_cooling=PS2020
	with_subgrid_chemistry=EAGLE
	with_subgrid_tracers=EAGLE
	with_subgrid_entropy_floor=EAGLE
	with_subgrid_stars=EAGLE
	with_subgrid_star_formation=EAGLE
	with_subgrid_feedback=EAGLE
	with_subgrid_black_holes=EAGLE
	with_subgrid_sink=none
	with_subgrid_extra_io=none
	enable_fof=yes
   ;;
   SPIN_JET_EAGLE)
	with_subgrid_cooling=EAGLE
	with_subgrid_chemistry=EAGLE
	with_subgrid_tracers=EAGLE
	with_subgrid_entropy_floor=EAGLE
	with_subgrid_stars=EAGLE
	with_subgrid_star_formation=EAGLE
	with_subgrid_feedback=EAGLE
	with_subgrid_black_holes=SPIN_JET
	with_subgrid_sink=none
	with_subgrid_extra_io=none
	enable_fof=yes
   ;;
   *)
      AC_MSG_ERROR([Unknown subgrid choice: $with_subgrid])
   ;;
esac

# Check if FoF is on.
AC_ARG_ENABLE([fof],
   [AS_HELP_STRING([--enable-fof],
     [Activate the friends-of-friends (FoF) code.],
   )],
   [enable_fof="$enableval"],
   [enable_fof="no"]
)
if test "$enable_fof" = "yes"; then
   AC_DEFINE([WITH_FOF], 1, [Enable FoF])
fi

# Check if stand-alone FoF is on.
AC_ARG_ENABLE([stand-alone-fof],
   [AS_HELP_STRING([--enable-stand-alone-fof],
     [Activate the compilation of the stand-alone friends-of-friends (FoF) post-processing tool.],
   )],
   [enable_standalone_fof="$enableval"],
   [enable_standalone_fof="no"]
)
if test "$enable_standalone_fof" = "yes"; then
   enable_fof="yes + stand-alone tool"
   AC_DEFINE([WITH_FOF], 1, [Enable FoF])
   AC_DEFINE([WITH_STAND_ALONE_FOF], 1, [Enable stand-alone FoF])
fi
AM_CONDITIONAL([HAVESTANDALONEFOF],[test $enable_standalone_fof = "yes"])

# Gravity scheme.
AC_ARG_WITH([gravity],
   [AS_HELP_STRING([--with-gravity=<scheme>],
      [Gravity scheme to use @<:@basic, with-multi-softening default: with-multi-softening@:>@]
   )],
   [with_gravity="$withval"],
   [with_gravity="with-multi-softening"]
)

case "$with_gravity" in
   with-potential)
      AC_MSG_ERROR([The gravity 'with-potential' scheme does not exist anymore. Please use the basic scheme which now contains potentials.])
   ;;
   with-multi-softening)
      AC_DEFINE([MULTI_SOFTENING_GRAVITY], [1], [Gravity scheme with per-particle type softening value and background particles])
   ;;
   basic)
      AC_DEFINE([DEFAULT_GRAVITY], [1], [Basic gravity scheme])
   ;;
   *)
      AC_MSG_ERROR([Unknown gravity scheme: $with_gravity])
   ;;
esac

AC_ARG_ENABLE([gravitational-potential],
   [AS_HELP_STRING([--disable-gravitational-potential],
     [Disable calculation of the gravitational potential.]
   )],
   [enable_gravitational_potential="$enableval"],
   [enable_gravitational_potential="yes"]
)
if test "$enable_gravitational_potential" = "no"; then
   AC_DEFINE([SWIFT_GRAVITY_NO_POTENTIAL],1,[Disable calculation of the gravitational potential])
fi

# Hydro scheme.
AC_ARG_WITH([hydro],
   [AS_HELP_STRING([--with-hydro=<scheme>],
      [Hydro dynamics to use @<:@gadget2, minimal, pressure-entropy, pressure-energy, pressure-energy-monaghan, phantom, gizmo-mfv, gizmo-mfm, shadowfax, planetary, sphenix, gasoline, anarchy-pu default: sphenix@:>@]
   )],
   [with_hydro="$withval"],
   [with_hydro="sphenix"]
)

case "$with_hydro" in
   none)
      AC_DEFINE([NONE_SPH], [1], [No hydro])
   ;;
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
   phantom)
      AC_DEFINE([PHANTOM_SPH], [1], [Phantom SPH])
   ;;
   gizmo-mfv)
      AC_DEFINE([GIZMO_MFV_SPH], [1], [GIZMO MFV SPH])
      need_riemann_solver=yes
      hydro_does_mass_flux=yes
   ;;
   gizmo-mfm)
      AC_DEFINE([GIZMO_MFM_SPH], [1], [GIZMO MFM SPH])
      need_riemann_solver=yes
   ;;
   shadowfax)
      AC_DEFINE([SHADOWFAX_SPH], [1], [Shadowfax SPH])
      need_riemann_solver=yes
   ;;
   planetary)
      AC_DEFINE([PLANETARY_SPH], [1], [Planetary SPH])
   ;;
   sphenix)
      AC_DEFINE([SPHENIX_SPH], [1], [SPHENIX SPH])
   ;;
   gasoline)
      AC_DEFINE([GASOLINE_SPH], [1], [Gasoline SPH])
   ;;
   anarchy-du)
      AC_DEFINE([SPHENIX_SPH], [1], [SPHENIX SPH])
   ;;
   anarchy-pu)
      AC_DEFINE([ANARCHY_PU_SPH], [1], [ANARCHY (PU) SPH])
   ;;

   *)
      AC_MSG_ERROR([Unknown hydrodynamics scheme: $with_hydro])
   ;;
esac

# SPMHD scheme.
AC_ARG_WITH([spmhd],
   [AS_HELP_STRING([--with-spmhd=<scheme>],
      [Magneto Hydro Dynamics SPH scheme to use @<:@none, direct-induction, direct-induction-fede, vector-potential default:none@:>@]
   )],
   [with_spmhd="$withval"],
   [with_spmhd="none"]
)

case "$with_spmhd" in
   none)
      AC_DEFINE([NONE_MHD], [1], [No mhd])
   ;;
   *)
      AC_MSG_ERROR([Unknown magneto-hydrodynamics scheme: $with_spmhd])
   ;;
esac

if test "$with_hydro" = "gizmo-mfm" -a "$with_spmhd" != "none"; then
  AC_MSG_ERROR([Cannot use an SPMHD scheme alongside a gizmo hydro solver!"])
fi
if test "$with_hydro" = "gizmo-mfv" -a "$with_spmhd" != "none"; then
  AC_MSG_ERROR([Cannot use an SPMHD scheme alongside a gizmo hydro solver!"])
fi
if test "$with_hydro" = "shadowfax" -a "$with_spmhd" != "none"; then
  AC_MSG_ERROR([Cannot use an SPMHD scheme alongside a gizmo hydro solver!"])
fi

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
    if test "$enable_debug_interactions_stars" = "yes"; then
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
    if test "$enable_debug_interactions" = "yes"; then
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

# Check if debugging interactions sinks is switched on.
AC_ARG_ENABLE([debug-interactions-sinks],
   [AS_HELP_STRING([--enable-debug-interactions-sinks],
     [Activate interaction debugging for sinks, logging a maximum of @<:@N@:>@ neighbours. Defaults to 256 if no value set.]
   )],
   [enable_debug_interactions_sinks="$enableval"],
   [enable_debug_interactions_sinks="no"]
)
if test "$enable_debug_interactions_sinks" != "no"; then
    AC_DEFINE([DEBUG_INTERACTIONS_SINKS],1,[Enable interaction debugging for sinks])
    if test "$enable_debug_interactions_sinks" = "yes"; then
      AC_DEFINE([MAX_NUM_OF_NEIGHBOURS_SINKS],256,[The maximum number of particle neighbours to be logged for sinks])
      [enable_debug_interactions_sinks="yes (Logging up to 256 neighbours)"]
    else
      AC_DEFINE_UNQUOTED([MAX_NUM_OF_NEIGHBOURS_SINKS], [$enableval] ,[The maximum number of particle neighbours to be logged for sinks])
      [enable_debug_interactions_sinks="yes (Logging up to $enableval neighbours)"]
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
      [equation of state @<:@ideal-gas, isothermal-gas, barotropic-gas, planetary default: ideal-gas@:>@]
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
     barotropic-gas)
       AC_DEFINE([EOS_BAROTROPIC_GAS], [1], [Barotropic gas equation of state])
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

# Adaptive softening
AC_ARG_WITH([adaptive-softening],
   [AS_HELP_STRING([--with-adaptive-softening=<yes/no>],
      [Adaptive softening @<:@no, yes, default: no@:>@]
   )],
   [with_adaptive_softening="$withval"],
   [with_adaptive_softening="no"]
)
case "$with_adaptive_softening" in
   no)
      AC_DEFINE([FIXED_SOFTENING], [1], [No adaptive softening])
   ;;
   yes)
      AC_DEFINE([ADAPTIVE_SOFTENING], [1], [Adaptive softening])
   ;;
   *)
      AC_MSG_ERROR([Unknown adaptive softening: $with_adaptive_softening])
   ;;
esac
# Verify that the configuration is allowed
if test "x$with_adaptive_softening" = "xyes" -a "$with_kernel" != "wendland-C2"; then
  AC_MSG_ERROR([Adaptive softening scheme requires the usage of the Wendland-C2 kernel!])
fi
if test "x$with_adaptive_softening" = "xyes" -a "$with_gravity" != "with-multi-softening"; then
  AC_MSG_ERROR([Adaptive softening scheme requires the usage of the multi-softening gravity scheme!])
fi


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

if test "x$need_riemann_solver" = "xyes" -a "$with_riemann" = "none"; then
  AC_MSG_ERROR([Hydro scheme $with_hydro requires selection of a Riemann solver!])
fi

# Hydro does mass flux?
if test "x$hydro_does_mass_flux" = "xyes"; then
    AC_DEFINE([HYDRO_DOES_MASS_FLUX], [1], [Hydro scheme with mass fluxes])
fi

#  chemistry function
AC_ARG_WITH([chemistry],
   [AS_HELP_STRING([--with-chemistry=<function>],
      [chemistry function @<:@none, GEAR_*, AGORA, QLA, EAGLE default: none@:>@
      For GEAR, you need to provide the number of elements (e.g. GEAR_10)]
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

with_chemistry_name="none"
case "$with_chemistry" in
   none)
      AC_DEFINE([CHEMISTRY_NONE], [1], [No chemistry function])
   ;;
   GEAR_*)
      AC_DEFINE([CHEMISTRY_GEAR], [1], [Chemistry taken from the GEAR model])
      number_element=${with_chemistry#*_}
      AC_DEFINE_UNQUOTED([GEAR_CHEMISTRY_ELEMENT_COUNT], [$number_element], [Number of element to follow])
      with_chemistry_name="GEAR (with $number_element elements)"
      with_chemistry="GEAR"
   ;;
   AGORA)
      AC_DEFINE([CHEMISTRY_AGORA], [1], [Chemistry taken from the AGORA model])
      with_chemistry_name="AGORA"
   ;;
   QLA)
      AC_DEFINE([CHEMISTRY_QLA], [1], [Chemistry taken from the Quick-Lyman-alpha model])
      with_chemistry_name="QLA"
   ;;
   EAGLE)
      AC_DEFINE([CHEMISTRY_EAGLE], [1], [Chemistry taken from the EAGLE model])
      with_chemistry_name="EAGLE (9 elements + smoothing)"
   ;;
   *)
      AC_MSG_ERROR([Unknown chemistry function: $with_chemistry])
   ;;
esac

if test "$with_chemistry" != "none"; then
   if test "$enable_hand_vec" = "yes"; then
      if test "$enable_vec" = "yes"; then
         if test "$with_hydro" = "gadget2"; then
            AC_MSG_ERROR([Cannot run with hand vectorisation and chemistry yet. Please use --disable-hand-vec])
         fi
      fi
   fi
fi

#  Cooling function
AC_ARG_WITH([cooling],
   [AS_HELP_STRING([--with-cooling=<model>],
      [cooling function @<:@none, const-du, const-lambda, QLA, EAGLE, PS2020, Wiersma, grackle_* default: none@:>@.
      For Grackle, you need to provide the primordial chemistry parameter (e.g. grackle_0)]
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

with_cooling_name="none"
case "$with_cooling" in
   none)
      AC_DEFINE([COOLING_NONE], [1], [No cooling function])
   ;;
   const-du)
      AC_DEFINE([COOLING_CONST_DU], [1], [Const du/dt cooling function])
      with_cooling_name=$with_cooling
   ;;
   const-lambda)
      AC_DEFINE([COOLING_CONST_LAMBDA], [1], [Const Lambda cooling function])
      with_cooling_name=$with_cooling
   ;;
   grackle_*)

      if test "$have_grackle" != "yes"; then
        AC_MSG_ERROR([Grackle cooling: You need the grackle library for Grackle cooling. (--with-grackle=PATH)])
      fi

      AC_DEFINE([COOLING_GRACKLE], [1], [Cooling via the grackle library])
      primordial_chemistry=${with_cooling#*_}
      AC_DEFINE_UNQUOTED([COOLING_GRACKLE_MODE], [$primordial_chemistry], [Grackle chemistry network])
      with_cooling_name="Grackle $primordial_chemistry"
      with_cooling="grackle"
   ;;
   QLA)
      AC_DEFINE([COOLING_QLA], [1], [Cooling following the Quick-Lyman-alpha model])
      with_cooling_name="QLA (Ploeckinger+20 tables) with constant primordial Z"
   ;;
   QLA-EAGLE)
      AC_DEFINE([COOLING_QLA_EAGLE], [1], [Cooling following the Quick-Lyman-alpha model])
      with_cooling_name="QLA (Wiersma+09 tables) with constant primordial Z"
   ;;
   EAGLE)
      AC_DEFINE([COOLING_EAGLE], [1], [Cooling following the EAGLE model (Wiersma+09 tables)])
      with_cooling_name="EAGLE (Wiersma+09 tables)"
   ;;
   Wiersma)
      AC_DEFINE([COOLING_EAGLE], [1], [Cooling following the EAGLE model (Wiersma+09 tables)])
      with_cooling_name="EAGLE (Wiersma+09 tables)"
   ;;
   PS2020)
      AC_DEFINE([COOLING_PS2020], [1], [Cooling following the PS2020 model (Ploeckinger+20 tables)])
      with_cooling_name="PS2020 (Ploeckinger+20 tables)"
   ;;
   COLIBRE)
      AC_MSG_ERROR([This cooling model changed name. Please use 'PS2020'.])
   ;;
   *)
      AC_MSG_ERROR([Unknown cooling function: $with_cooling])
   ;;
esac

if test "$with_cooling" = "EAGLE" || test "$with_cooling" = "PS2020"; then
   if test "$with_chemistry" = "none"; then
      AC_MSG_ERROR([Cannot run with EAGLE or PS2020 cooling without chemistry. Please pick a chemistry model])
   fi
fi

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

#  Extra fields added to snapshots at i/o time
AC_ARG_WITH([extra_io],
   [AS_HELP_STRING([--with-extra-io=<function>],
      [Extra i/o field: @<:@none, EAGLE default: none@:>@]
   )],
   [with_extra_io="$withval"],
   [with_extra_io="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_extra_io" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-extra_io together])
   else
      with_extra_io="$with_subgrid_extra_io"
   fi
fi

case "$with_extra_io" in
   none)
      AC_DEFINE([EXTRA_IO_NONE], [1], [No extra_io function])
   ;;
   EAGLE)
      AC_DEFINE([EXTRA_IO_EAGLE], [1], [Extra i/o fields taken from the EAGLE model])
   ;;
   *)
      AC_MSG_ERROR([Unknown extra-io choice: $with_extra_io])
   ;;
esac

# Stellar model.
AC_ARG_WITH([stars],
   [AS_HELP_STRING([--with-stars=<model>],
      [Stellar model to use @<:@none, basic, EAGLE, GEAR, default: basic@:>@]
   )],
   [with_stars="$withval"],
   [with_stars="basic"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_stars" != "basic"; then
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
   basic)
      AC_DEFINE([STARS_BASIC], [1], [Basic stellar model])
   ;;
   none)
      AC_DEFINE([STARS_NONE], [1], [No stellar model])
   ;;
   *)
      AC_MSG_ERROR([Unknown stellar model: $with_stars])
   ;;
esac

# Feedback model
AC_ARG_WITH([feedback],
   [AS_HELP_STRING([--with-feedback=<model>],
      [Feedback model to use @<:@none, EAGLE, EAGLE-thermal, EAGLE-kinetic, GEAR, AGORA default: none@:>@]
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

with_feedback_name="none"
case "$with_feedback" in
   EAGLE-kinetic)
      AC_DEFINE([FEEDBACK_EAGLE_KINETIC], [1], [EAGLE kinetic stellar feedback and evolution model])
      with_feedback_name="EAGLE kinetic stellar feedback and evolution model"
   ;;
   EAGLE-thermal)
      AC_DEFINE([FEEDBACK_EAGLE_THERMAL], [1], [EAGLE thermal stellar feedback and evolution model])
      with_feedback_name="EAGLE thermal stellar feedback and evolution model"
   ;;
   EAGLE)
      AC_DEFINE([FEEDBACK_EAGLE_THERMAL], [1], [EAGLE thermal stellar feedback and evolution model])
      with_feedback_name="EAGLE thermal stellar feedback and evolution model"
   ;;
   GEAR)
      AC_DEFINE([FEEDBACK_GEAR], [1], [GEAR stellar feedback and evolution model])
      with_feedback_name="GEAR"
   ;;
   AGORA)
      AC_DEFINE([FEEDBACK_AGORA], [1], [AGORA stellar feedback and evolution model])
      with_feedback_name="AGORA"
   ;;
   none)
      AC_DEFINE([FEEDBACK_NONE], [1], [No feedback])
   ;;

   *)
      AC_MSG_ERROR([Unknown feedback model: $with_feedback])
   ;;
esac

# For models using ray-based feedback, let the users set different numbers of rays
AC_ARG_WITH([number-of-SNII-rays],
	    [AS_HELP_STRING([--with-number-of-SNII-rays], [Number of rays to use for the SNII feedback (default: 1)])],
	    [with_number_of_SNII_rays="$withval"],
	    [with_number_of_SNII_rays="1"])
AC_DEFINE_UNQUOTED([FEEDBACK_NR_RAYS_SNII], [$with_number_of_SNII_rays], [Number of rays to use for the SNII feedback])

AC_ARG_WITH([number-of-SNIa-rays],
	    [AS_HELP_STRING([--with-number-of-SNIa-rays], [Number of rays to use for the SNIa feedback (default: 1)])],
	    [with_number_of_SNIa_rays="$withval"],
	    [with_number_of_SNIa_rays="1"])
AC_DEFINE_UNQUOTED([FEEDBACK_NR_RAYS_SNIa], [$with_number_of_SNIa_rays], [Number of rays to use for the SNIa feedback])

AC_ARG_WITH([number-of-AGN-rays],
	    [AS_HELP_STRING([--with-number-of-AGN-rays], [Number of rays to use for the AGN feedback (default: 50)])],
	    [with_number_of_AGN_rays="$withval"],
	    [with_number_of_AGN_rays="50"])
AC_DEFINE_UNQUOTED([FEEDBACK_NR_RAYS_AGN], [$with_number_of_AGN_rays], [Number of rays to use for the AGN feedback])


# Black hole model.
AC_ARG_WITH([black-holes],
   [AS_HELP_STRING([--with-black-holes=<model>],
      [Black holes model to use @<:@none, EAGLE, SPIN_JET default: none@:>@]
   )],
   [with_black_holes="$withval"],
   [with_black_holes="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_black_holes" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-black-holes together])
   else
      with_black_holes="$with_subgrid_black_holes"
   fi
fi

case "$with_black_holes" in
   none)
      AC_DEFINE([BLACK_HOLES_NONE], [1], [No black hole model])
   ;;
   EAGLE)
      AC_DEFINE([BLACK_HOLES_EAGLE], [1], [EAGLE black hole model])
   ;;
   SPIN_JET)
      AC_DEFINE([BLACK_HOLES_SPIN_JET], [1], [Spin and jet black hole model])
      with_black_holes="SPIN_JETS (Husko+22)"
   ;;
   *)
      AC_MSG_ERROR([Unknown black-hole model: $with_black_holes])
   ;;
esac

# Sink model.
AC_ARG_WITH([sink],
   [AS_HELP_STRING([--with-sink=<model>],
      [Sink particle model to use @<:@none, GEAR, default: none@:>@]
   )],
   [with_sink="$withval"],
   [with_sink="none"]
)

if test "$with_subgrid" != "none"; then
   if test "$with_sink" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-sink together])
   else
      with_sink="$with_subgrid_sink"
   fi
fi

case "$with_sink" in
   none)
      AC_DEFINE([SINK_NONE], [1], [No sink particle model])
   ;;
   GEAR)
    AC_DEFINE([SINK_GEAR], [1], [GEAR sink particle model])
    ;;
    *)
      AC_MSG_ERROR([Unknown sink particle model model: $with_sink])
   ;;
esac

# Forcing terms
AC_ARG_WITH([forcing],
   [AS_HELP_STRING([--with-forcing=<term>],
      [Hydrodynamics forcing terms @<:@none, roberts-flow, roberts-flow-acceleration default: none@:>@]
   )],
   [with_forcing="$withval"],
   [with_forcing="none"]
)
case "$with_forcing" in
   none)
      AC_DEFINE([FORCING_NONE], [1], [No external forcing terms])
   ;;
   roberts-flow)
      AC_DEFINE([FORCING_ROBERTS_FLOW], [1], [Roberts' flow external forcing terms])
   ;;
   roberts-flow-acceleration)
      AC_DEFINE([FORCING_ROBERTS_FLOW_ACCELERATION], [1], [Roberts' flow external forcing terms entering the equations as an acceleration term])
   ;;
   *)
      AC_MSG_ERROR([Unknown external forcing term: $with_forcing])
   ;;
esac


#  External potential
AC_ARG_WITH([ext-potential],
   [AS_HELP_STRING([--with-ext-potential=<pot>],
      [external potential @<:@none, point-mass, point-mass-softened, isothermal, nfw, nfw-mn, hernquist, hernquist-sdmh05, disc-patch, sine-wave, MWPotential2014, constant, default: none@:>@]
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
   hernquist-sdmh05)
      AC_DEFINE([EXTERNAL_POTENTIAL_HERNQUIST_SDMH05], [1], [Hernquist external potential following Springel, Di Matteo & Hernquist 2005])
   ;;
   nfw)
      AC_DEFINE([EXTERNAL_POTENTIAL_NFW], [1], [Navarro-Frenk-White external potential])
   ;;
   nfw-mn)
      AC_DEFINE([EXTERNAL_POTENTIAL_NFW_MN], [1], [Navarro-Frenk-White + Miyamoto-Nagai disk external potential])
   ;;
   disc-patch)
      AC_DEFINE([EXTERNAL_POTENTIAL_DISC_PATCH], [1], [Disc-patch external potential])
   ;;
   sine-wave)
      AC_DEFINE([EXTERNAL_POTENTIAL_SINE_WAVE], [1], [Sine wave external potential in 1D])
   ;;
   point-mass-softened)
      AC_DEFINE([EXTERNAL_POTENTIAL_POINTMASS_SOFT], [1], [Softened point-mass potential with form 1/(r^2 + softening^2).])
   ;;
   MWPotential2014)
      AC_DEFINE([EXTERNAL_POTENTIAL_MWPotential2014], [1], [Milky-Way like potential composed of a Navarro-Frenk-White + Miyamoto-Nagai disk + Power spherical cuttoff external potential.])
   ;;
   constant)
      AC_DEFINE([EXTERNAL_POTENTIAL_CONSTANT], [1], [Constant gravitational acceleration.])
   ;;
   *)
      AC_MSG_ERROR([Unknown external potential: $with_potential])
   ;;
esac

#  Entropy floor
AC_ARG_WITH([entropy-floor],
    [AS_HELP_STRING([--with-entropy-floor=<floor>],
       [entropy floor @<:@none, QLA, EAGLE, default: none@:>@]
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
   QLA)
      AC_DEFINE([ENTROPY_FLOOR_QLA], [1], [Quick Lyman-alpha entropy floor])
   ;;
   EAGLE)
      AC_DEFINE([ENTROPY_FLOOR_EAGLE], [1], [EAGLE entropy floor])
   ;;
   *)
      AC_MSG_ERROR([Unknown entropy floor model])
   ;;
esac

#  Pressure floor
AC_ARG_WITH([pressure-floor],
    [AS_HELP_STRING([--with-pressure-floor=<floor>],
       [pressure floor @<:@none, GEAR, default: none@:>@
       The hydro model needs to be compatible.]
    )],
    [with_pressure_floor="$withval"],
    [with_pressure_floor="none"]
)
if test "$with_subgrid" != "none"; then
   if test "$with_pressure_floor" != "none"; then
      AC_MSG_ERROR([Cannot provide with-subgrid and with-pressure-floor together])
   else
      with_pressure_floor="$with_subgrid_pressure_floor"
   fi
fi

case "$with_pressure_floor" in
   none)
      AC_DEFINE([PRESSURE_FLOOR_NONE], [1], [No pressure floor])
   ;;
   GEAR)
      AC_DEFINE([PRESSURE_FLOOR_GEAR], [1], [GEAR pressure floor])
   ;;
   *)
      AC_MSG_ERROR([Unknown pressure floor model])
   ;;
esac

#  Star formation
AC_ARG_WITH([star-formation],
    [AS_HELP_STRING([--with-star-formation=<sfm>],
       [star formation @<:@none, QLA, EAGLE, GEAR, default: none@:>@]
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
   QLA)
      AC_DEFINE([STAR_FORMATION_QLA], [1], [Quick Lyman-alpha star formation model)])
   ;;
   EAGLE)
      AC_DEFINE([STAR_FORMATION_EAGLE], [1], [EAGLE star formation model (Schaye and Dalla Vecchia (2008))])
   ;;
   GEAR)
      AC_DEFINE([STAR_FORMATION_GEAR], [1], [GEAR star formation model (Revaz and Jablonka (2018))])
   ;;
   *)
      AC_MSG_ERROR([Unknown star formation model])
   ;;
esac

AC_ARG_WITH([gadget2-physical-constants],
    [AS_HELP_STRING([--with-gadget2-physical-constants],
       [Use the same physical constants (G, Msun, Mpc) as the Gadget-2 code rather than more up-to-date values.]
    )],
    [with_gadget2_physical_constants="$withval"],
    [with_gadget2_physical_constants="no"]
)
if test "$with_gadget2_physical_constants" = "yes"; then
   AC_DEFINE([SWIFT_USE_GADGET2_PHYSICAL_CONSTANTS],1,[Use the same physical constants as Gadget-2])
fi

#  Gravity multipole order
AC_ARG_WITH([multipole-order],
   [AS_HELP_STRING([--with-multipole-order=<order>],
      [order of the multipole and gravitational field expansion @<:@ default: 4@:>@]
   )],
   [with_multipole_order="$withval"],
   [with_multipole_order="4"]
)
AC_DEFINE_UNQUOTED([SELF_GRAVITY_MULTIPOLE_ORDER], [$with_multipole_order], [Multipole order])

#  Radiative transfer scheme
AC_ARG_WITH([rt],
   [AS_HELP_STRING([--with-rt=<scheme>],
      [Radiative transfer scheme to use @<:@none, GEAR_*, SPHM1RT_*, debug default: none@:>@.
      For GEAR and SPHM1RT, the number of photon groups (e.g. GEAR_4) needs to be provided.]
   )],
   [with_rt="$withval"],
   [with_rt="none"]
)

# For GEAR-RT scheme: Select a RT Riemann solver
AC_ARG_WITH([rt-riemann-solver],
   [AS_HELP_STRING([--with-rt-riemann-solver=<scheme>],
      [Riemann solver for the moments of the ratiadiative transfer equation with the M1 closure to use @<:@none, HLL, GLF, default: none@:>@.
      For the GEAR RT scheme, you need to select one Riemann solver.]
   )],
   [with_rt_riemann_solver="$withval"],
   [with_rt_riemann_solver="none"]
)

case "$with_rt_riemann_solver" in
    none)
        AC_DEFINE([RT_RIEMANN_SOLVER_NONE], [1], [No RT riemann solver])
        ;;
    HLL)
        AC_DEFINE([RT_RIEMANN_SOLVER_HLL], [1], [HLL RT riemann solver])
        ;;
    GLF)
        AC_DEFINE([RT_RIEMANN_SOLVER_GLF], [1], [GLF RT riemann solver])
        ;;
    *)
        AC_MSG_ERROR("Unknown RT Riemann solver: $with_rt_riemann_solver")
        ;;
esac

case "$with_rt" in
   none)
      AC_DEFINE([RT_NONE], [1], [No radiative transfer scheme])
   ;;
   GEAR_*)
      AC_DEFINE([RT_GEAR], [1], [GEAR M1 closure scheme])
      number_group=${with_rt#*_}
      AC_DEFINE_UNQUOTED([RT_NGROUPS], [$number_group], [Number of photon groups to follow])

      if test "$number_group" = "0"; then
          AC_MSG_ERROR([GEAR-RT: Cannot work with zero photon groups])
      fi
      if ! test $number_group -eq $number_group; then
          # abuse -eq to check whether $number_group is an integer. -eq
          # only works with those.
          AC_MSG_ERROR([GEAR-RT: Cannot work with non-integer photon groups])
      fi

      if test "$enable_debugging_checks" = "yes"; then
          AC_DEFINE([SWIFT_RT_DEBUG_CHECKS], [1], [additional debugging checks for RT])
      fi

      if test "$with_hydro" != "gizmo-mfv"; then
          AC_MSG_ERROR([GEAR-RT: Cannot work without gizmo-mfv hydro. Compile using --with-hydro=gizmo-mfv])
      fi

      if test "$with_rt_riemann_solver" = "none"; then
          AC_MSG_ERROR([GEAR-RT: You need to select an RT Riemann solver (--with-rt-riemann-solver=...)])
      fi

      if test "$have_grackle" != "yes"; then
        AC_MSG_ERROR([GEAR-RT: You need the grackle library for GEAR-RT. (--with-grackle=PATH)])
      fi

   ;;
   debug)
      AC_DEFINE([RT_DEBUG], [1], [debugging scheme])
      AC_DEFINE([SWIFT_RT_DEBUG_CHECKS], [1], [additional debugging checks for RT])
   ;;
   SPHM1RT_*)
      AC_DEFINE([RT_SPHM1RT], [1], [SPHM1RT radiative transfer scheme (Chan+21: 2102.08404)])
      number_group=${with_rt#*_}
      AC_DEFINE_UNQUOTED([RT_NGROUPS], [$number_group], [Number of photon groups to follow])

      if test "$number_group" -eq "0"; then
          AC_MSG_ERROR([SPHM1RT: Cannot work with zero photon groups])
      fi
      if test $number_group -lt 4; then
          AC_MSG_ERROR([SPHM1RT: Cannot work with < four photon groups for now])
      fi
      if ! test $number_group -eq $number_group; then
          # abuse -eq to check whether $number_group is an integer. -eq
          # only works with those.
          AC_MSG_ERROR([SPHM1RT: Cannot work with non-integer photon groups])
      fi
      AC_MSG_CHECKING([for Sundials libraries])
      AC_MSG_RESULT($have_sundials)
      if test "$have_sundials" != "yes"; then
         AC_MSG_ERROR([The Sundials library is not present. Sundials is required for the SPHM1RT module.])
      fi
   ;;
   *)
      AC_MSG_ERROR([Unknown radiative transfer scheme: $with_rt])
   ;;
esac


# Check for git, needed for revision stamps.
AC_PATH_PROG([GIT_CMD], [git])
AC_SUBST([GIT_CMD])

# Make the documentation. Add conditional to handle disable option.
DX_DOXYGEN_FEATURE(OFF)
DX_INIT_DOXYGEN(SWIFT, doc/Doxyfile, doc/)
AM_CONDITIONAL([HAVE_DOXYGEN], [test "$ac_cv_path_ac_pt_DX_DOXYGEN" != ""])

# Check if using EAGLE extra I/O
AM_CONDITIONAL([HAVEEAGLEEXTRAIO], [test "${with_extra_io}" = "EAGLE"])

# Check if using QLA cooling
AM_CONDITIONAL([HAVEQLACOOLING], [test "$with_cooling" = "QLA"])
AM_CONDITIONAL([HAVEQLAEAGLECOOLING], [test "$with_cooling" = "QLA-EAGLE"])

# Check if using EAGLE cooling
AM_CONDITIONAL([HAVEEAGLECOOLING], [test "$with_cooling" = "EAGLE"])

# Check if using PS2020 cooling
AM_CONDITIONAL([HAVEPS2020COOLING], [test "$with_cooling" = "PS2020"])

# Check if using EAGLE feedback
AM_CONDITIONAL([HAVEEAGLETHERMALFEEDBACK], [test "${with_feedback%-thermal}" = "EAGLE"])
AM_CONDITIONAL([HAVEEAGLEKINETICFEEDBACK], [test "$with_feedback" = "EAGLE-kinetic"])

# check if using grackle cooling
AM_CONDITIONAL([HAVEGRACKLECOOLING], [test "$with_cooling" = "grackle"])

# check if using gear feedback
AM_CONDITIONAL([HAVEGEARFEEDBACK], [test "$with_feedback" = "GEAR"])

# check if using gear feedback
AM_CONDITIONAL([HAVEAGORAFEEDBACK], [test "$with_feedback" = "AGORA"])

# check if using SPHENIX
AM_CONDITIONAL([HAVE_SPHENIX], [test "$with_hydro" = "sphenix"])

# check if using GADGET2 SPH
AM_CONDITIONAL([HAVE_GADGET2], [test "$with_hydro" = "gadget2"])

# check if using none chemistry
AM_CONDITIONAL([HAVE_CHEMISTRY_NONE], [test "$with_chemistry" = "none"])

# check if using GEAR chemistry
AM_CONDITIONAL([HAVE_CHEMISTRY_GEAR], [test "$with_chemistry" = "GEAR" || test "$with_chemistry" = "GEAR_DIFFUSION"])

# check if using AGORA chemistry
AM_CONDITIONAL([HAVE_CHEMISTRY_AGORA], [test "$with_chemistry" = "AGORA" || test "$with_chemistry" = "GEAR_DIFFUSION"])

# check if using default stars
AM_CONDITIONAL([HAVE_STARS_BASIC], [test "$with_stars" = "basic"])

# check if using GEAR stars
AM_CONDITIONAL([HAVE_STARS_GEAR], [test "$with_stars" = "GEAR"])

# check if using default star formation
AM_CONDITIONAL([HAVE_STAR_FORMATION_DEFAULT], [test "$with_star_formation" = "none"])

# check if using GEAR star formation
AM_CONDITIONAL([HAVE_STAR_FORMATION_GEAR], [test "$with_star_formation" = "GEAR"])

# check if using multi softening gravity
AM_CONDITIONAL([HAVE_GRAVITY_MULTISOFTENING], [test "$with_gravity" = "with-multi-softening"])

# Check if using SPHM1RT radiative transfer
AM_CONDITIONAL([HAVESPHM1RTRT], [test "${with_rt:0:7}" = "SPHM1RT"])

# Check if using GEAR-RT radiative transfer
AM_CONDITIONAL([HAVEGEARRT], [test "${with_rt:0:4}" = "GEAR"])




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
AC_CONFIG_FILES([tests/testNeutrinoCosmology.sh], [chmod +x tests/testNeutrinoCosmology.sh])
AC_CONFIG_FILES([tests/output_list_params.yml])

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
    - threaded/openmp   : $have_threaded_fftw / $have_openmp_fftw
    - MPI               : $have_mpi_fftw
    - ARM               : $have_arm_fftw
   GSL enabled          : $have_gsl
   HEALPix C enabled    : $have_chealpix
   libNUMA enabled      : $have_numa
   GRACKLE enabled      : $have_grackle
   Sundials enabled     : $have_sundials
   Special allocators   : $have_special_allocator
   CPU profiler         : $have_profiler
   Pthread barriers     : $have_pthread_barrier
   VELOCIraptor enabled : $have_velociraptor
   FoF activated:       : $enable_fof
   Lightcones enabled   : $enable_lightcone

   Hydro scheme       : $with_hydro
   Dimensionality     : $with_dimension
   Kernel function    : $with_kernel
   Equation of state  : $with_eos
   Adiabatic index    : $with_gamma
   Riemann solver     : $with_riemann
   SPMHD scheme       : $with_spmhd
   Adaptive softening : $with_adaptive_softening

   Gravity scheme      : $with_gravity
   Multipole order     : $with_multipole_order
   Compute potential   : $enable_gravitational_potential
   No gravity below ID : $no_gravity_below_id
   Make gravity glass  : $gravity_glass_making
   External potential  : $with_potential
   Forcing terms       : $with_forcing

   Pressure floor       : $with_pressure_floor
   Entropy floor        : $with_entropy_floor
   Cooling function     : $with_cooling_name
   Chemistry            : $with_chemistry_name
   Tracers              : $with_tracers
   Stellar model        : $with_stars
   Star formation model : $with_star_formation
   Star feedback model  : $with_feedback_name
   Sink particle model  : $with_sink
   Black holes model    : $with_black_holes
   Radiative transfer   : $with_rt
   Extra i/o            : $with_extra_io

   Atomic operations in tasks  : $enable_atomics_within_tasks
   Individual timers           : $enable_timers
   Task debugging              : $enable_task_debugging
   Threadpool debugging        : $enable_threadpool_debugging
   Debugging checks            : $enable_debugging_checks
   Interaction debugging       : $enable_debug_interactions
   Stars interaction debugging : $enable_debug_interactions_stars
   Naive interactions          : $enable_naive_interactions
   Naive stars interactions    : $enable_naive_interactions_stars
   Gravity checks              : $gravity_force_checks
   Custom icbrtf               : $enable_custom_icbrtf
   Boundary particles          : $boundary_particles
   Fixed boundary particles    : $fixed_boundary_particles
   Planetary fixed entropy     : $planetary_fixed_entropy
   Ghost statistics            : $ghost_stats

   Continuous Sim. Data Stream : $with_csds

 ------------------------]
)

#  Report any unrecognised options loudly.
if test "$ac_unrecognized_opts" != ""; then
AC_MSG_RESULT([
  *** WARNING: unrecognized options: $ac_unrecognized_opts ***
])
fi
