# ===========================================================================
#        http://www.gnu.org/software/autoconf-archive/ax_lib_hdf5.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_HDF5([serial/parallel])
#
# DESCRIPTION
#
#   This macro provides tests of the availability of HDF5 library.
#
#   The optional macro argument should be either 'serial' or 'parallel'. The
#   former only looks for serial HDF5 installations via h5cc. The latter
#   only looks for parallel HDF5 installations via h5pcc. If the optional
#   argument is omitted, serial installations will be preferred over
#   parallel ones.
#
#   The macro adds a --with-hdf5 option accepting one of three values:
#
#     no   - do not check for the HDF5 library.
#     yes  - do check for HDF5 library in standard locations.
#     path - complete path to where lib/libhdf5* libraries and
#            include/H5* include files reside.
#
#   If HDF5 is successfully found, this macro calls
#
#     AC_SUBST(HDF5_VERSION)
#     AC_SUBST(HDF5_CFLAGS)
#     AC_SUBST(HDF5_CPPFLAGS)
#     AC_SUBST(HDF5_LDFLAGS)
#     AC_SUBST(HDF5_FFLAGS)
#     AC_SUBST(HDF5_FLIBS)
#     AC_DEFINE(HAVE_HDF5)
#
#   and sets with_hdf5="yes".  Additionally, the macro sets
#   with_hdf5_fortran="yes" if a matching Fortran wrapper script is found.
#   Note that Autconf's Fortran support is not used to perform this check.
#   H5CC and H5FC will contain the appropriate serial or parallel HDF5
#   wrapper script locations.
#
#   If HDF5 is disabled or not found, this macros sets with_hdf5="no" and
#   with_hdf5_fortran="no".
#
#   Your configuration script can test $with_hdf to take any further
#   actions. HDF5_{C,CPP,LD}FLAGS may be used when building with C or C++.
#   HDF5_F{FLAGS,LIBS} should be used when building Fortran applications.
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for HDF5 support
#        AX_LIB_HDF5()
#
#     2) dnl Check for serial HDF5 support
#        AX_LIB_HDF5([serial])
#
#     3) dnl Check for parallel HDF5 support
#        AX_LIB_HDF5([parallel])
#
#   One could test $with_hdf5 for the outcome or display it as follows
#
#     echo "HDF5 support:  $with_hdf5"
#
# LICENSE
#
#   Copyright (c) 2009 Timothy Brown <tbrown@freeshell.org>
#   Copyright (c) 2010 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 6

AC_DEFUN([AX_LIB_HDF5], [

AC_REQUIRE([AC_PROG_SED])
AC_REQUIRE([AC_PROG_AWK])
AC_REQUIRE([AC_PROG_GREP])

dnl Check first argument is one of the recognized values.
dnl Fail eagerly if is incorrect as this simplifies case statements below.
if   test "m4_normalize(m4_default([$1],[]))" = ""        ; then
    : # Recognized value
elif test "m4_normalize(m4_default([$1],[]))" = "serial"  ; then
    : # Recognized value
elif test "m4_normalize(m4_default([$1],[]))" = "parallel"; then
    : # Recognized value
else
    AC_MSG_ERROR([
Unrecognized value for AX[]_LIB_HDF5 within configure.ac.
If supplied, argument 1 must be either 'serial' or 'parallel'.
])
fi

dnl Add a default --with-hdf5 configuration option.
AC_ARG_WITH([hdf5],
  AS_HELP_STRING(
    [--with-hdf5=[yes/no/PATH]],
    m4_case(m4_normalize([$1]),
            [serial],   [location of h5cc for serial HDF5 configuration],
            [parallel], [location of h5pcc for parallel HDF5 configuration],
            [location of h5cc or h5pcc for HDF5 configuration])
  ),
  [if test "$withval" = "no"; then
     with_hdf5="no"
   elif test "$withval" = "yes"; then
     with_hdf5="yes"
   else
     with_hdf5="yes"
     H5CC="$withval"
   fi],
   [with_hdf5="yes"]
)

dnl Set defaults to blank
HDF5_VERSION=""
HDF5_CFLAGS=""
HDF5_CPPFLAGS=""
HDF5_LDFLAGS=""
HDF5_FFLAGS=""
HDF5_FLIBS=""

dnl Try and find hdf5 compiler tools and options.
if test "$with_hdf5" = "yes"; then
    if test -z "$H5CC"; then
        dnl Check to see if H5CC is in the path.
        AC_PATH_PROGS(
            [H5CC],
            m4_case(m4_normalize([$1]),
                [serial],   [h5cc],
                [parallel], [h5pcc],
                [h5cc h5pcc]),
            [])
    else
        AC_MSG_CHECKING([Using provided HDF5 C wrapper])
        AC_MSG_RESULT([$H5CC])
    fi
    AC_MSG_CHECKING([for HDF5 libraries])
    if test ! -x "$H5CC"; then
        AC_MSG_RESULT([no])
        AC_MSG_WARN(m4_case(m4_normalize([$1]),
            [serial],  [
Unable to locate serial HDF5 compilation helper script 'h5cc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5cc.
HDF5 support is being disabled (equivalent to --with-hdf5=no).
],            [parallel],[
Unable to locate parallel HDF5 compilation helper script 'h5pcc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5pcc.
HDF5 support is being disabled (equivalent to --with-hdf5=no).
],            [
Unable to locate HDF5 compilation helper scripts 'h5cc' or 'h5pcc'.
Please specify --with-hdf5=<LOCATION> as the full path to h5cc or h5pcc.
HDF5 support is being disabled (equivalent to --with-hdf5=no).
]))
        with_hdf5="no"
        with_hdf5_fortran="no"
    else
        dnl h5cc provides both AM_ and non-AM_ options
        dnl depending on how it was compiled either one of
        dnl these are empty. Lets roll them both into one.

        dnl Look for "HDF5 Version: X.Y.Z"
        HDF5_VERSION=$(eval $H5CC -showconfig | grep 'HDF5 Version:' \
            | $AWK '{print $[]3}')

dnl A ideal situation would be where everything we needed was
dnl in the AM_* variables. However most systems are not like this
dnl and seem to have the values in the non-AM variables.
dnl
dnl We try the following to find the flags:
dnl (1) Look for "NAME:" tags
dnl (2) Look for "NAME/H5_NAME:" tags
dnl (3) Look for "AM_NAME:" tags
dnl
        dnl (1)
        dnl Look for "CFLAGS: "
        HDF5_CFLAGS=$(eval $H5CC -showconfig | grep '\bCFLAGS:'   \
            | $AWK -F: '{print $[]2}')
        dnl Look for "CPPFLAGS"
        HDF5_CPPFLAGS=$(eval $H5CC -showconfig | grep '\bCPPFLAGS:' \
            | $AWK -F: '{print $[]2}')
        dnl Look for "FFLAGS: "
        HDF5_FFLAGS=$(eval $H5CC -showconfig | grep '\bFFLAGS:'   \
            | $AWK -F: '{print $[]2}')
        dnl Look for "LD_FLAGS"
        HDF5_LDFLAGS=$(eval $H5CC -showconfig | grep '\bLDFLAGS:'   \
            | $AWK -F: '{print $[]2}')

        dnl (2)
        dnl CFLAGS/H5_CFLAGS: .../....
        dnl We could use $SED with something like the following
        dnl 's/CFLAGS.*\/H5_CFLAGS.*[:]\(.*\)\/\(.*\)/\1/p'
        if test -z "$HDF5_CFLAGS"; then
            HDF5_CFLAGS=$(eval $H5CC -showconfig \
                | $SED -n 's/CFLAGS.*[:]\(.*\)\/\(.*\)/\1/p')
        fi
        dnl Look for "CPPFLAGS"
        if test -z "$HDF5_CPPFLAGS"; then
            HDF5_CPPFLAGS=$(eval $H5CC -showconfig \
                | $SED -n 's/CPPFLAGS.*[:]\(.*\)\/\(.*\)/\1/p')
        fi
        dnl Look for "FFLAGS"
        if test -z "$HDF5_FFLAGS"; then
            HDF5_FFLAGS=$(eval $H5CC -showconfig \
                | $SED -n 's/FFLAGS.*[:]\(.*\)\/\(.*\)/\1/p')
        fi
        dnl Look for "LD_FLAGS"
        if test -z "$HDF5_LDFLAGS"; then
            HDF5_LDFLAGS=$(eval $H5CC -showconfig \
                | $SED -n 's/LDFLAGS.*[:]\(.*\)\/\(.*\)/\1/p')
        fi

        dnl (3)
        dnl Check to see if these are not empty strings. If so
        dnl find the AM_ versions and use them.
        if test -z "$HDF5_CFLAGS"; then
            HDF5_CFLAGS=$(eval $H5CC -showconfig \
                | grep 'AM_CFLAGS:' | $AWK -F: '{print $[]2}')
        fi
        if test -z "$HDF5_CPPFLAGS"; then
            HDF5_CPPFLAGS=$(eval $H5CC -showconfig \
                | grep 'AM_CPPFLAGS:' | $AWK -F: '{print $[]2}')
        fi
        dnl No equivalent AM_FFLAGS test in h5cc -showconfig
        if test -z "$HDF5_LDFLAGS"; then
            HDF5_LDFLAGS=$(eval $H5CC -showconfig \
                | grep 'AM_LDFLAGS:' | $AWK -F: '{print $[]2}')
        fi

        dnl Frustratingly, the necessary -Idir,-Ldir still may not be found!
        dnl Attempt to pry any more required include directories from wrapper.
        for arg in `$H5CC -c -show` # Wrapper compile mode
        do
          case "$arg" in #(
            -I*) echo $HDF5_CPPFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_CPPFLAGS="$arg $HDF5_CPPFLAGS"
              ;;
          esac
        done
        for arg in `$H5CC -show` # Wrapper link mode
        do
          case "$arg" in #(
            -L*) echo $HDF5_LDFLAGS | $GREP -e "$arg" 2>&1 >/dev/null \
                  || HDF5_LDFLAGS="$arg $HDF5_LDFLAGS"
              ;;
          esac
        done

        AC_MSG_RESULT([yes (version $[HDF5_VERSION])])

        dnl Look for any extra libraries also needed to link properly
        EXTRA_LIBS=$(eval $H5CC -showconfig | grep 'Extra libraries:'\
            | $AWK -F: '{print $[]2}')

        dnl Look for HDF5's high level library
        ax_lib_hdf5_save_LDFLAGS=$LDFLAGS
        ax_lib_hdf5_save_LIBS=$LIBS
        LDFLAGS=$HDF5_LDFLAGS
        AC_HAVE_LIBRARY([hdf5_hl],
                        [HDF5_LDFLAGS="$HDF5_LDFLAGS -lhdf5_hl"],
                        [],
                        [-lhdf5 $EXTRA_LIBS])
        LIBS=$ax_lib_hdf5_save_LIBS
        LDFLAGS=$ax_lib_hdf5_save_LDFLAGS

        dnl Add the HDF5 library itself
        HDF5_LDFLAGS="$HDF5_LDFLAGS -lhdf5"

        dnl Add any EXTRA_LIBS afterwards
        if test "$EXTRA_LIBS"; then
            HDF5_LDFLAGS="$HDF5_LDFLAGS $EXTRA_LIBS"
        fi

        AC_MSG_CHECKING([for matching HDF5 Fortran wrapper])
        dnl Presume HDF5 Fortran wrapper is just a name variant from H5CC
        H5FC=$(eval echo -n $H5CC | $SED -n 's/cc$/fc/p')
        if test -x "$H5FC"; then
            AC_MSG_RESULT([$H5FC])
            with_hdf5_fortran="yes"
            AC_SUBST([H5FC])

            dnl Again, pry any remaining -Idir/-Ldir from compiler wrapper
            for arg in `$H5FC -show`
            do
              case "$arg" in #(
                -I*) echo $HDF5_FFLAGS | $GREP -e "$arg" >/dev/null \
                      || HDF5_FFLAGS="$arg $HDF5_FFLAGS"
                  ;;#(
                -L*) echo $HDF5_FFLAGS | $GREP -e "$arg" >/dev/null \
                      || HDF5_FFLAGS="$arg $HDF5_FFLAGS"
                     dnl HDF5 installs .mod files in with libraries,
                     dnl but some compilers need to find them with -I
                     echo $HDF5_FFLAGS | $GREP -e "-I${arg#-L}" >/dev/null \
                      || HDF5_FFLAGS="-I${arg#-L} $HDF5_FFLAGS"
                  ;;
              esac
            done

            dnl Make Fortran link line by inserting Fortran libraries
            for arg in $HDF5_LDFLAGS
            do
              case "$arg" in #(
                -lhdf5_hl) HDF5_FLIBS="$HDF5_FLIBS -lhdf5hl_fortran $arg"
                  ;; #(
                -lhdf5)    HDF5_FLIBS="$HDF5_FLIBS -lhdf5_fortran $arg"
                  ;; #(
                *) HDF5_FLIBS="$HDF5_FLIBS $arg"
                  ;;
              esac
            done
        else
            AC_MSG_RESULT([no])
            with_hdf5_fortran="no"
        fi

	AC_SUBST([HDF5_VERSION])
	AC_SUBST([HDF5_CFLAGS])
	AC_SUBST([HDF5_CPPFLAGS])
	AC_SUBST([HDF5_LDFLAGS])
	AC_SUBST([HDF5_FFLAGS])
	AC_SUBST([HDF5_FLIBS])
	AC_DEFINE([HAVE_HDF5], [1], [Defined if you have HDF5 support])
    fi
fi

])
