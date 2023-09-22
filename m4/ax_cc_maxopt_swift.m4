# ===========================================================================
#       https://www.gnu.org/software/autoconf-archive/ax_cc_maxopt.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CC_MAXOPT_SWIFT
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
#   SWIFT version does not add the optimization flags to CFLAGS.  This allows
#   us to use these flags in automake so we have have different levels of
#   optimization on a library by library basis. It and also means the user
#   CFLAGS can just extended these rather than needing to copy them all.
#
#   The flags are returned in the variable OPT_CFLAGS. This is set regardless
#   of the CFLAGS value, unlike the original version of this macro.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2008 Matteo Frigo
#   Copyright (c) 2023 Peter W. Draper
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

#serial 18

AC_DEFUN([AX_CC_MAXOPT_SWIFT],
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AX_COMPILER_VENDOR])
AC_REQUIRE([AC_CANONICAL_HOST])

AC_ARG_ENABLE(portable-binary, [AS_HELP_STRING([--enable-portable-binary], [disable compiler optimizations that would produce unportable binaries])],
	acx_maxopt_portable=$enableval, acx_maxopt_portable=no)

# Try to determine "good" native compiler flags.
case $ax_cv_c_compiler_vendor in
  dec) OPT_CFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host"
       if test "x$acx_maxopt_portable" = xno; then
          OPT_CFLAGS="$OPT_CFLAGS -arch host"
       fi
     ;;

  sun) OPT_CFLAGS="-native -fast -xO5 -dalign"
       if test "x$acx_maxopt_portable" = xyes; then
	  CFLAGS="$OPT_CFLAGS -xarch=generic"
       fi
     ;;

  hp)  CFLAGS="+Oall +Optrs_ansi +DSnative"
       if test "x$acx_maxopt_portable" = xyes; then
          CFLAGS="$OPT_CFLAGS +DAportable"
       fi
     ;;

  ibm) if test "x$acx_maxopt_portable" = xno; then
          xlc_opt="-qarch=auto -qtune=auto"
       else
          xlc_opt="-qtune=auto"
       fi
       AX_CHECK_COMPILE_FLAG($xlc_opt,
                OPT_CFLAGS="-O3 -qansialias -w $xlc_opt",
               [OPT_CFLAGS="-O3 -qansialias -w"
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

  intel) OPT_CFLAGS="-O3 -ansi-alias"
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
                    *0?6[[78ab]]?:*:*:*|?6[[78ab]]?:*:*:*|6[[78ab]]?:*:*:*) icc_flags="-xK" ;;
                    *0?6[[9d]]?:*:*:*|?6[[9d]]?:*:*:*|6[[9d]]?:*:*:*|*1?65?:*:*:*) icc_flags="-xSSE2 -xB -xK" ;;
                    *0?6e?:*:*:*|?6e?:*:*:*|6e?:*:*:*) icc_flags="-xSSE3 -xP -xO -xB -xK" ;;
                    *0?6f?:*:*:*|?6f?:*:*:*|6f?:*:*:*|*1?66?:*:*:*) icc_flags="-xSSSE3 -xT -xB -xK" ;;
                    *1?6[[7d]]?:*:*:*) icc_flags="-xSSE4.1 -xS -xT -xB -xK" ;;
                    *1?6[[aef]]?:*:*:*|*2?6[[5cef]]?:*:*:*) icc_flags="-xSSE4.2 -xS -xT -xB -xK" ;;
                    *2?6[[ad]]?:*:*:*) icc_flags="-xAVX -SSE4.2 -xS -xT -xB -xK" ;; # Sandy-bridge
                    *3?6[[ae]]?:*:*:*) icc_flags="-xCORE-AVX-I -xAVX -SSE4.2 -xS -xT -xB -xK" ;; #Ivy-bridge
                    *3?6[[cf]]?:*:*:*|*4?6[[56]]?:*:*:*|*4?6[[ef]]?:*:*:*) icc_flags="-xCORE-AVX2 -xCORE-AVX-I -xAVX -SSE4.2 -xS -xT -xB -xK" ;; # Haswell
                    *3?6d?:*:*:*|*4?6[[7f]]?:*:*:*|*5?66?:*:*:*) icc_flags="-xCORE-AVX2 -xCORE-AVX-I -xAVX -SSE4.2 -xS -xT -xB -xK" ;; # Broadwell
                    *4?6[[de]]?:*:*:*) icc_flags="-xCORE-AVX2 -xCORE-AVX-I -xAVX -SSE4.2 -xS -xT -xB -xK" ;; # Skylake
                    *5?6[[56]]?:*:*:*) icc_flags="-xCORE-AVX512 -xCORE-AVX2 -xCORE-AVX-I -xAVX -SSE4.2 -xS -xT -xB -xK" ;; # Skylake-AVX512
                    *5?67?:*:*:*) icc_flags="-xMIC-AVX512 -xCORE-AVX2 -xCORE-AVX-I -xAVX -SSE4.2 -xS -xT -xB -xK" ;; # Knights-Landing
                    *8?6[[de]]?:*:*:*|*9?6[[de]]?:*:*:*) icc_flags="-xCORE-AVX2 -xCORE-AVX-I -xAVX -SSE4.2 -xS -xT -xB -xK" ;;# Kabylake
                    *000?f[[346]]?:*:*:*|?f[[346]]?:*:*:*|f[[346]]?:*:*:*) icc_flags="-xSSE3 -xP -xO -xN -xW -xK" ;;
                    *00??f??:*:*:*|??f??:*:*:*|?f??:*:*:*|f??:*:*:*) icc_flags="-xSSE2 -xN -xW -xK" ;;
                  esac ;;
                *:68747541:444d4163:69746e65) # AMDs with AVX2 support.
                  case $ax_cv_gcc_x86_cpuid_1 in
                    *061?f??:*:*:*|61?f??:*:*:*) icc_flags="-march=core-avx2" ;;
                    *06??f??:*:*:*|6??f??:*:*:*) icc_flags="-march=core-avx2" ;;
                    *070?f??:*:*:*|70?f??:*:*:*) icc_flags="-march=core-avx2" ;;
                                   83?f??:*:*:*) icc_flags="-march=core-avx2"
                                                 OPT_CFLAGS="$OPT_CFLAGS -fma -ftz -fomit-frame-pointer";; # EPYC ROME
                                   a0?f??:*:*:*) icc_flags="-march=core-avx2"
                                                 OPT_CFLAGS="$OPT_CFLAGS -fma -ftz -fomit-frame-pointer";; # MILAN

                  esac ;;
              esac ;;
          esac
          if test "x$icc_flags" != x; then
            for flag in $icc_flags; do
              AX_CHECK_COMPILE_FLAG($flag, [icc_archflag=$flag; break])
            done
          fi
          AC_MSG_CHECKING([for icc architecture flag])
          AC_MSG_RESULT($icc_archflag)
          if test "x$icc_archflag" != xunknown; then
            OPT_CFLAGS="$OPT_CFLAGS $icc_archflag"
          fi
        fi
     ;;

  clang)
     # default optimization flags for clang on all systems
     OPT_CFLAGS="-O3 -fomit-frame-pointer"

     # Always good optimisation to have
     AX_CHECK_COMPILE_FLAG(-fstrict-aliasing, [OPT_CFLAGS="$OPT_CFLAGS -fstrict-aliasing"])

     # note that we enable "unsafe" fp optimization with other compilers, too
     AX_CHECK_COMPILE_FLAG(-ffast-math, [OPT_CFLAGS="$OPT_CFLAGS -ffast-math"])

     # not all codes will benefit from this.
     AX_CHECK_COMPILE_FLAG(-funroll-loops, [OPT_CFLAGS="$OPT_CFLAGS -funroll-loops"])

     AX_GCC_ARCHFLAG($acx_maxopt_portable, [OPT_CFLAGS="$OPT_CFLAGS $ax_cv_gcc_archflag"])
     ;;

  gnu)
     # default optimization flags for gcc on all systems
     OPT_CFLAGS="-O3 -fomit-frame-pointer"

     # -malign-double for x86 systems
     AX_CHECK_COMPILE_FLAG(-malign-double, [OPT_CFLAGS="$OPT_CFLAGS -malign-double"])

     #  -fstrict-aliasing for gcc-2.95+
     AX_CHECK_COMPILE_FLAG(-fstrict-aliasing, [OPT_CFLAGS="$OPT_CFLAGS -fstrict-aliasing"])

     # note that we enable "unsafe" fp optimization with other compilers, too
     AX_CHECK_COMPILE_FLAG(-ffast-math, [OPT_CFLAGS="$OPT_CFLAGS -ffast-math"])

     # not all codes will benefit from this.
     AX_CHECK_COMPILE_FLAG(-funroll-loops, [OPT_CFLAGS="$OPT_CFLAGS -funroll-loops"])

     AX_GCC_ARCHFLAG($acx_maxopt_portable, [OPT_CFLAGS="$OPT_CFLAGS $ax_cv_gcc_archflag"])
     ;;

  microsoft)
     # default optimization flags for MSVC opt builds
     OPT_CFLAGS="-O2"
     ;;
esac

if test -z "$OPT_CFLAGS"; then
   echo ""
   echo "********************************************************"
   echo "* WARNING: Don't know the best optimizing flags for    *"
   echo "* this system a default of -O3 will be used.           *"
   echo "* Use ./configure --disable-optimization CFLAGS=...    *"
   echo "* to add your own flags.                               *"
   echo "********************************************************"
   echo ""
   OPT_CFLAGS="-O3"
fi

AX_CHECK_COMPILE_FLAG($OPT_CFLAGS, [], [
   echo ""
   echo "********************************************************************"
   echo "* WARNING: The guessed optimization flags don't seem to work with  *"
   echo "* your compiler.                                       *"
   echo "* Use ./configure CFLAGS=... to specify your own flags *"
   echo "********************************************************"
   echo ""
   OPT_CFLAGS=""
])

])
