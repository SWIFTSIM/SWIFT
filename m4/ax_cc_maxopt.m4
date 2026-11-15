# ===========================================================================
#       https://www.gnu.org/software/autoconf-archive/ax_cc_maxopt.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CC_MAXOPT
#
# DESCRIPTION
#
#   Try to turn on "good" C optimization flags for various compilers and
#   architectures, for some definition of "good". (In our case, good for
#   FFTW and hopefully for other scientific codes. Modify as needed.)
#
#   The user can override the flags by setting the CFLAGS environment
#   variable. The user can also specify --enable-portable-binary in order to
#   disable any optimization flags that might result in a binary that only
#   runs on the host architecture.
#
#   Note also that the flags assume that ANSI C aliasing rules are followed
#   by the code (e.g. for gcc's -fstrict-aliasing), and that floating-point
#   computations can be re-ordered as needed.
#
#   Requires macros: AX_CHECK_COMPILE_FLAG, AX_COMPILER_VENDOR,
#   AX_GCC_ARCHFLAG, AX_GCC_X86_CPUID.
#
#   NOTE: this copy has diverged from the Autoconf Archive and should not be
#   replaced wholesale by a newer upstream revision without merging what
#   follows. Measured against upstream serial 23, it adds branches for clang,
#   for icx (oneAPI), and for nvc, Cray and Fujitsu, none of which upstream
#   selects useful flags for; an Intel -x table that reaches current CPUs
#   rather than stopping at Haswell, with AMD parts given -march= instead
#   because the -x codes gate on a GenuineIntel check at run time; and a
#   preference for asking the compiler to target the build host itself, with
#   -march=native, -mcpu=native, -xHost or -fast, leaving the CPUID tables as
#   the fallback for portable and cross builds. It also spells the Intel
#   aliasing flag -ansi-alias rather than -ansi_alias, and tests
#   ac_test_CFLAGS against "set" rather than against the empty string, which
#   is the form configure.ac in this tree depends on.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 18 (modified for SWIFT)

AC_DEFUN([AX_CC_MAXOPT],
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AX_COMPILER_VENDOR])
AC_REQUIRE([AC_CANONICAL_HOST])

AC_ARG_ENABLE(portable-binary, [AS_HELP_STRING([--enable-portable-binary], [disable compiler optimizations that would produce unportable binaries])],
	acx_maxopt_portable=$enableval, acx_maxopt_portable=no)

# Try to determine "good" native compiler flags if none specified via CFLAGS
if test "$ac_test_CFLAGS" != "set"; then
  case $ax_cv_c_compiler_vendor in
    dec) CFLAGS="$CFLAGS -newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host"
	 if test "x$acx_maxopt_portable" = xno; then
           CFLAGS="$CFLAGS -arch host"
         fi;;

    sun) CFLAGS="$CFLAGS -native -fast -xO5 -dalign"
	 if test "x$acx_maxopt_portable" = xyes; then
	   CFLAGS="$CFLAGS -xarch=generic"
         fi;;

    hp)  CFLAGS="$CFLAGS +Oall +Optrs_ansi +DSnative"
	 if test "x$acx_maxopt_portable" = xyes; then
	   CFLAGS="$CFLAGS +DAportable"
	 fi;;

    ibm) if test "x$acx_maxopt_portable" = xno; then
           xlc_opt="-qarch=auto -qtune=auto"
	 else
           xlc_opt="-qtune=auto"
	 fi
         AX_CHECK_COMPILE_FLAG($xlc_opt,
		CFLAGS="$CFLAGS -O3 -qansialias -w $xlc_opt",
               [CFLAGS="$CFLAGS -O3 -qansialias -w"
                echo "******************************************************"
                echo "*  You seem to have the IBM  C compiler.  It is      *"
                echo "*  recommended for best performance that you use:    *"
                echo "*                                                    *"
                echo "*    CFLAGS=-O3 -qarch=xxx -qtune=xxx -qansialias -w *"
                echo "*                      ^^^        ^^^                *"
                echo "*  where xxx is pwr2, pwr3, 604, or whatever kind of *"
                echo "*  CPU you have.  (Set the CFLAGS environment var.   *"
                echo "*  and re-run configure.)  For more info, man cc.    *"
                echo "******************************************************"])
         ;;

    intel | oneapi) CFLAGS="$CFLAGS -O3 -ansi-alias"
	if test "x$acx_maxopt_portable" = xno; then
	  icc_archflag=unknown
	  icc_flags=""
	  case $host_cpu in
	    i686*|x86_64*)
              # icc accepts gcc assembly syntax, so these should work:
	      AX_GCC_X86_CPUID(0)
              AX_GCC_X86_CPUID(1)
	      case $ax_cv_gcc_x86_cpuid_0 in # see AX_GCC_ARCHFLAG
                *:756e6547:6c65746e:49656e69) # Intel
                  case $ax_cv_gcc_x86_cpuid_1 in
		    # The table below this comment stops at Kaby Lake and Skylake-AVX512, so
		    # every more recent CPU fell through it with icc_flags empty and got no
		    # -x flag at all, leaving both compilers on their SSE2 default. Unlike the
		    # gcc path there is no AX_EXT rescue: configure.ac skips SIMD_FLAGS for the
		    # Intel vendor. These entries are grouped by instruction set rather than by
		    # microarchitecture, since that is what -x actually selects.
		    # Atom cores come first: they are SSE4.2-only, and model 0x4d collides with
		    # the Skylake pattern further down, which would hand an Avoton -xCORE-AVX2.
		    *3?67?:*:*:*|*[[45]]?6[[acd]]?:*:*:*|*5?6[[cf]]?:*:*:*|*7?6[[5a]]?:*:*:*|*8?6[[6a]]?:*:*:*|*9?6[[6c]]?:*:*:*) icc_flags="-xSSE4.2" ;; # Silvermont..Tremont
		    *5?65?:*:*:*|*6?6[[6ac]]?:*:*:*|*7?6[[de]]?:*:*:*|*9?6d?:*:*:*|*8?6[[cdf]]?:*:*:*|*a?6[[7de]]?:*:*:*|*c?6f?:*:*:*|*4??f??:*:*:*) icc_flags="-xCORE-AVX512 -xCORE-AVX2 -xCORE-AVX-I -xAVX -xSSE4.2" ;; # AVX-512 parts
		    *3?6[[cdf]]?:*:*:*|*4?6[[567ef]]?:*:*:*|*5?6[[6e]]?:*:*:*|*8?6e?:*:*:*|*9?6[[7ae]]?:*:*:*|*a?6[[56acf]]?:*:*:*|*b?6[[567adef]]?:*:*:*|*c?6[[56c]]?:*:*:*|*d?6[[57d]]?:*:*:*) icc_flags="-xCORE-AVX2 -xCORE-AVX-I -xAVX -xSSE4.2" ;; # Alder Lake onwards
		    *8?65?:*:*:*) icc_flags="-xMIC-AVX512 -xCORE-AVX2 -xAVX -xSSE4.2" ;; # Knights Mill
		    *0?6[[78ab]]?:*:*:*|?6[[78ab]]?:*:*:*|6[[78ab]]?:*:*:*) icc_flags="-xK" ;;
		    *0?6[[9d]]?:*:*:*|?6[[9d]]?:*:*:*|6[[9d]]?:*:*:*|*1?65?:*:*:*) icc_flags="-xSSE2 -xB -xK" ;;
		    *0?6e?:*:*:*|?6e?:*:*:*|6e?:*:*:*) icc_flags="-xSSE3 -xP -xO -xB -xK" ;;
		    *0?6f?:*:*:*|?6f?:*:*:*|6f?:*:*:*|*1?66?:*:*:*) icc_flags="-xSSSE3 -xT -xB -xK" ;;
		    *1?6[[7d]]?:*:*:*) icc_flags="-xSSE4.1 -xS -xT -xB -xK" ;;
		    *1?6[[aef]]?:*:*:*|*2?6[[5cef]]?:*:*:*) icc_flags="-xSSE4.2 -xS -xT -xB -xK" ;;
		    *2?6[[ad]]?:*:*:*) icc_flags="-xAVX -xSSE4.2 -xS -xT -xB -xK" ;; # Sandy-bridge
		    *3?6[[ae]]?:*:*:*) icc_flags="-xCORE-AVX-I -xAVX -xSSE4.2 -xS -xT -xB -xK" ;; #Ivy-bridge
		    *3?6[[cf]]?:*:*:*|*4?6[[56]]?:*:*:*|*4?6[[ef]]?:*:*:*) icc_flags="-xCORE-AVX2 -xCORE-AVX-I -xAVX -xSSE4.2 -xS -xT -xB -xK" ;; # Haswell
		    *3?6d?:*:*:*|*4?6[[7f]]?:*:*:*|*5?66?:*:*:*) icc_flags=" -xCORE-AVX2 -xCORE-AVX-I -xAVX -xSSE4.2 -xS -xT -xB -xK" ;; # Broadwell
		    *4?6[[de]]?:*:*:*) icc_flags="-xCORE-AVX2 -xCORE-AVX-I -xAVX -xSSE4.2 -xS -xT -xB -xK" ;; # Skylake
		    *5?6[[56]]?:*:*:*) icc_flags="-xCORE-AVX512 -xCORE-AVX2 -xCORE-AVX-I -xAVX -xSSE4.2 -xS -xT -xB -xK" ;; # Skylake-AVX512
		    *5?67?:*:*:*) icc_flags="-xMIC-AVX512 -xCORE-AVX2 -xCORE-AVX-I -xAVX -xSSE4.2 -xS -xT -xB -xK" ;; # Knights-Landing
		    *8?6[[de]]?:*:*:*|*9?6[[de]]?:*:*:*) icc_flags="-xCORE-AVX2 -xCORE-AVX-I -xAVX -xSSE4.2 -xS -xT -xB -xK" ;;# Kabylake
		    *000?f[[346]]?:*:*:*|?f[[346]]?:*:*:*|f[[346]]?:*:*:*) icc_flags="-xSSE3 -xP -xO -xN -xW -xK" ;;
		    *00??f??:*:*:*|??f??:*:*:*|?f??:*:*:*|f??:*:*:*) icc_flags="-xSSE2 -xN -xW -xK" ;;
		    # Unknown recent Intel: extended model >= 9 postdates Skylake, so AVX2 is a
		    # safe floor. The AVX2-less Tremont parts carry extended model 8 and 9 but
		    # are matched explicitly above, so they cannot reach this. Without it an
		    # unrecognised CPU gets no -x flag at all and falls back to SSE2.
		    *[[9a-f]]?6??:*:*:*) icc_flags="-xCORE-AVX2 -xCORE-AVX-I -xAVX -xSSE4.2" ;;
                  esac ;;
                *:68747541:444d4163:69746e65) # AMDs with AVX2 support.
                  case $ax_cv_gcc_x86_cpuid_1 in
                    *061?f??:*:*:*|61?f??:*:*:*) icc_flags="-march=core-avx2" ;;
                    *06??f??:*:*:*|6??f??:*:*:*) icc_flags="-march=core-avx2" ;;
                    *070?f??:*:*:*|70?f??:*:*:*) icc_flags="-march=core-avx2" ;;
                                   83?f??:*:*:*) icc_flags="-march=core-avx2"
                                                 CFLAGS="$CFLAGS -fma -ftz -fomit-frame-pointer";; # ROME
                                   a0?f??:*:*:*) icc_flags="-march=core-avx2"
                                                 CFLAGS="$CFLAGS -fma -ftz -fomit-frame-pointer";; # MILAN
                                   a1?f??:*:*:*) icc_flags="-axCORE-AVX512"
                                                 CFLAGS="$CFLAGS -march=skylake-avx512 -fma -ftz -fomit-frame-pointer";; # GENOA
                                   aa?f??:*:*:*) icc_flags="-axCORE-AVX512"
                                                 CFLAGS="$CFLAGS -march=skylake-avx512 -fma -ftz -fomit-frame-pointer";; # BERGAMO

                  esac ;;
              esac ;;
          esac
          # icx (oneAPI) dropped the single-letter and pre-SSE4.2 processor
          # codes that classic icc still accepts, so do not offer it flags its
          # driver will reject. Any host new enough to be worth building with
          # icx reaches one of the instruction-set entries above.
          if test "x$ax_cv_c_compiler_vendor" = xoneapi; then
            icc_oneapi_flags=""
            for flag in $icc_flags; do
              case $flag in
                -xSSE4.2|-xAVX|-xCORE-AVX-I|-xCORE-AVX2|-xCORE-AVX512|-xMIC-AVX512)
                  icc_oneapi_flags="$icc_oneapi_flags $flag" ;;
              esac
            done
            icc_flags=$icc_oneapi_flags
          fi
          # -xHost targets the build host directly, which is more accurate and
          # stays current without this table being edited, so offer it ahead of
          # the table flags; the loop below falls through to them if the driver
          # rejects it. Intel parts only: on AMD the -x codes gate on a
          # GenuineIntel check at run time, which is why the AMD entries above
          # use -march= instead. When cross compiling the cpuid is unknown, so
          # this does not match and -xHost is not offered.
          case $ax_cv_gcc_x86_cpuid_0 in
            *:756e6547:6c65746e:49656e69) icc_flags="-xHost $icc_flags" ;;
          esac
          if test "x$icc_flags" != x; then
            for flag in $icc_flags; do
              AX_CHECK_COMPILE_FLAG($flag, [icc_archflag=$flag; break])
            done
          fi
          AC_MSG_CHECKING([for icc architecture flag])
	  AC_MSG_RESULT($icc_archflag)
          if test "x$icc_archflag" != xunknown; then
            CFLAGS="$CFLAGS $icc_archflag"
          fi
        fi
	;;

    clang)
     # default optimization flags for clang on all systems
     CFLAGS="$CFLAGS -O3 -fomit-frame-pointer"

     # Always good optimisation to have
     AX_CHECK_COMPILE_FLAG(-fstrict-aliasing, CFLAGS="$CFLAGS -fstrict-aliasing")

     # note that we enable "unsafe" fp optimization with other compilers, too
     AX_CHECK_COMPILE_FLAG(-ffast-math, CFLAGS="$CFLAGS -ffast-math")

     # not all codes will benefit from this.
     AX_CHECK_COMPILE_FLAG(-funroll-loops, CFLAGS="$CFLAGS -funroll-loops")

     # Prefer the compiler's own host detection to the CPUID tables in
     # AX_GCC_ARCHFLAG. It is more accurate, it stays current without this file
     # being edited, and on Darwin it is the only thing that can work at all.
     # Skip it when a portable binary is asked for, and when cross compiling,
     # where native would describe the build machine and not the target. The
     # tables remain the fallback for both of those cases, as does an explicit
     # --with-gcc-arch=<arch>, which must keep overriding the guess.
     ax_maxopt_gotnative=no
     if test "x$acx_maxopt_portable" = xno && test "x$cross_compiling" = xno \
        && test -z "$with_gcc_arch"; then
       case $host_cpu in
         aarch64*|arm64*|powerpc*) ax_maxopt_nativeflag="-mcpu=native" ;;
         *) ax_maxopt_nativeflag="-march=native" ;;
       esac
       AX_CHECK_COMPILE_FLAG([$ax_maxopt_nativeflag],
         [CFLAGS="$CFLAGS $ax_maxopt_nativeflag"; ax_maxopt_gotnative=yes])
     fi
     if test "x$ax_maxopt_gotnative" = xno; then
       AX_GCC_ARCHFLAG($acx_maxopt_portable)
     fi
     ;;

    gnu)
     # default optimization flags for gcc on all systems
     CFLAGS="$CFLAGS -O3 -fomit-frame-pointer"

     # -malign-double for x86 systems
     AX_CHECK_COMPILE_FLAG(-malign-double, CFLAGS="$CFLAGS -malign-double")

     #  -fstrict-aliasing for gcc-2.95+
     AX_CHECK_COMPILE_FLAG(-fstrict-aliasing,
        CFLAGS="$CFLAGS -fstrict-aliasing")

     # note that we enable "unsafe" fp optimization with other compilers, too
     AX_CHECK_COMPILE_FLAG(-ffast-math, CFLAGS="$CFLAGS -ffast-math")

     # not all codes will benefit from this.
     AX_CHECK_COMPILE_FLAG(-funroll-loops, CFLAGS="$CFLAGS -funroll-loops")

     # Prefer the compiler's own host detection to the CPUID tables in
     # AX_GCC_ARCHFLAG. It is more accurate, it stays current without this file
     # being edited, and on Darwin it is the only thing that can work at all.
     # Skip it when a portable binary is asked for, and when cross compiling,
     # where native would describe the build machine and not the target. The
     # tables remain the fallback for both of those cases, as does an explicit
     # --with-gcc-arch=<arch>, which must keep overriding the guess.
     ax_maxopt_gotnative=no
     if test "x$acx_maxopt_portable" = xno && test "x$cross_compiling" = xno \
        && test -z "$with_gcc_arch"; then
       case $host_cpu in
         aarch64*|arm64*|powerpc*) ax_maxopt_nativeflag="-mcpu=native" ;;
         *) ax_maxopt_nativeflag="-march=native" ;;
       esac
       AX_CHECK_COMPILE_FLAG([$ax_maxopt_nativeflag],
         [CFLAGS="$CFLAGS $ax_maxopt_nativeflag"; ax_maxopt_gotnative=yes])
     fi
     if test "x$ax_maxopt_gotnative" = xno; then
       AX_GCC_ARCHFLAG($acx_maxopt_portable)
     fi
     ;;

    portland | nvhpc)
     # nvc, the NVIDIA HPC SDK compiler, formerly PGI. Both vendor strings are
     # matched on purpose: nvc defines __NVCOMPILER as well as __PGI, and while
     # the AX_COMPILER_VENDOR in this tree only tests for the latter, upstream
     # added an nvhpc entry ahead of portland. Accepting either means a refresh
     # of that macro cannot silently drop this branch. Its reference guide says
     # of -fast that "the appropriate -tp option is automatically included to
     # enable generation of code optimized for the type of system on which
     # compilation is performed", so the one flag covers both the optimisation
     # level and the host targeting that -march=native gives elsewhere. That
     # also makes it the wrong choice for a portable binary, hence the test.
     if test "x$acx_maxopt_portable" = xno && test "x$cross_compiling" = xno; then
       AX_CHECK_COMPILE_FLAG(-fast, CFLAGS="$CFLAGS -fast", [CFLAGS="$CFLAGS -O3"])
     else
       CFLAGS="$CFLAGS -O3"
     fi
     ;;

    cray)
     # Classic Cray C only. CCE 9 and later are clang based and define
     # __clang__, which AX_COMPILER_VENDOR tests before _CRAYC, so those are
     # handled by the clang branch above. No architecture flag is set here: on
     # a Cray the cc wrapper takes the target from the loaded craype-* module,
     # and overriding that from configure is more likely to fight it than help.
     CFLAGS="$CFLAGS -O3"
     ;;

    fujitsu)
     # fcc in Trad mode. In Clang mode it defines __clang__ and is handled by
     # the clang branch above. -Kfast is the aggregate optimisation flag and
     # already implies the relaxed floating point that the gcc and clang paths
     # ask for with -ffast-math, as well as targeting the build host.
     AX_CHECK_COMPILE_FLAG(-Kfast, CFLAGS="$CFLAGS -Kfast", [CFLAGS="$CFLAGS -O3"])
     ;;

    microsoft)
     # default optimization flags for MSVC opt builds
     CFLAGS="$CFLAGS -O2"
     ;;
  esac

  # Nothing above can work out the target when cross compiling. AX_GCC_ARCHFLAG
  # skips its tables, on every architecture and not just some, and asking the
  # compiler about a machine it is not running on is meaningless, so the native
  # paths are skipped too. The result is a build with no architecture flags at
  # all, which succeeds and is merely slower than it should be, so say so
  # rather than leave it to be discovered. --with-gcc-arch is exempt because it
  # replaces the detection instead of refining it, and so still applies here.
  if test "x$cross_compiling" = xyes && test -z "$with_gcc_arch"; then
     AC_MSG_WARN([cross compiling: no architecture flags have been selected, so
this build will not be tuned for the machine it is meant to run on. Use
--with-gcc-arch=<arch> with GCC or clang, or set CFLAGS yourself, to choose one.])
  fi

  if test -z "$CFLAGS"; then
	echo ""
	echo "********************************************************"
        echo "* WARNING: Don't know the best CFLAGS for this system  *"
        echo "* Use ./configure CFLAGS=... to specify your own flags *"
	echo "* (otherwise, a default of CFLAGS=-O3 will be used)    *"
	echo "********************************************************"
	echo ""
        CFLAGS="$CFLAGS -O3"
  fi

  AX_CHECK_COMPILE_FLAG($CFLAGS, [], [
	echo ""
        echo "********************************************************"
        echo "* WARNING: The guessed CFLAGS don't seem to work with  *"
        echo "* your compiler.                                       *"
        echo "* Use ./configure CFLAGS=... to specify your own flags *"
        echo "********************************************************"
        echo ""
  ])

fi
])
