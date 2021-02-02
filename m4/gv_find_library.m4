#
# SYNOPSIS
#
#   GV_FIND_LIBRARY(NAME, VARNAME, PKGNAME, LIBNAME, FUNC)
#
#     NAME    : Name to use in help strings
#     VARNAME : Defines macro HAVE_VARNAME and shell variable USE_VARNAME
#               if library found
#     PKGNAME : package name used by pkg-config
#     LIBNAME : Name used in the library file name
#     FUNC    : A function in the library
#
# DESCRIPTION
#
#   Attempts to set up the specified library using pkg-config
#
#   Uses pkg-config to try to find package PKGNAME. If successful defines
#   preprocessor macro HAVE_$VARNAME and sets USE_$VARNAME = "yes", otherwise
#   sets USE_$VARNAME = "no". 
#
#   Also sets ${VARNAME}_LIBS and ${VARNAME}_CFLAGS on success.
#
# LAST MODIFICATION
#
#   07/03/09 Version used in Gadgetviewer
#   02/02/21 Use export when temporarily setting PKG_CONFIG_PATH
#

AC_DEFUN([GV_FIND_LIBRARY],[

# Allow user to enable/disable library
AC_ARG_WITH([$1], AC_HELP_STRING([--with-$1@<:@=PATH@:>@],[use the $1 library]),
		  [USE_$2=$withval ; GV_SPEC="yes" ],[USE_$2="check" ; GV_SPEC="no" ])

# Figure out if we have been given a path
if test $USE_$2 = "yes" ; then
  GV_HAVE_$2_PATH="no"
  GV_$2_PATH=""
elif test $USE_$2 = "no" ; then
  GV_HAVE_$2_PATH="no"
  GV_$2_PATH=""
elif test $USE_$2 = "check" ; then
  GV_HAVE_$2_PATH="no"
  GV_$2_PATH=""
else
  GV_HAVE_$2_PATH="yes"
  GV_$2_PATH=$USE_$2
fi

GV_FOUND="no"

# Don't do anything if library has been disabled explicitly
if test $USE_$2 != "no" ; then

  # Add path to PKG_CONFIG_PATH if we have one
  TMP_PKG_CONFIG_PATH=$PKG_CONFIG_PATH
  if test $GV_HAVE_$2_PATH = "yes" ; then
    export PKG_CONFIG_PATH=$GV_$2_PATH/pkgconfig/:$GV_$2_PATH/lib/pkgconfig/:$PKG_CONFIG_PATH
  fi

  # Try to set it up with pkg-config
  GV_PKG="yes"
  PKG_CHECK_MODULES([$2], [$3], [], 
  			  [echo Unable to find $1 with pkg-config, will try to link to it anyway... ; GV_PKG="no"])

  # Restore original PKG_CONFIG_PATH
  export PKG_CONFIG_PATH=$TMP_PKG_CONFIG_PATH

  # If that didn't work and flags haven't been supplied by hand but we have a path, try sensible defaults
  if test ${GV_PKG} = "no" ; then
    if test ${GV_HAVE_$2_PATH} = "yes" ; then
      # CFLAGS
      if test X${$2_CFLAGS} = X ; then
        $2_CFLAGS=-I${GV_$2_PATH}/include/
      fi
      # LIBS
      if test X${$2_LIBS} = X ; then
        $2_LIBS="-L${GV_$2_PATH}/lib/"
      fi
    fi
  fi

  # Try to link to the library
  TMP_LDFLAGS=$LDFLAGS
  LDFLAGS="${LDFLAGS} ${$2_LIBS}"
  AC_CHECK_LIB([$4], [$5], [GV_FOUND="yes"], 
  		     [echo ; echo WARNING: Unable to link to $1 library. See config.log for details ; echo])
  LDFLAGS=$TMP_LDFLAGS

  if test $GV_FOUND = "no" ; then
    # If we can't link the library and it was explicitly asked for, abort	
    if test $GV_SPEC = "yes" ; then
      AC_MSG_ERROR([Unable to link to requested library: $1])
    fi
    # If the test failed, don't set flags
    $2_LIBS=""
    $2_CFLAGS=""
  else
    # If the test worked make sure -lwhatever is included if we didn't
    # use pkg-config
    if test $GV_PKG = "no" ; then
      $2_LIBS="${$2_LIBS} -l$4"
    fi
  fi

  AC_SUBST($2_LIBS)
  AC_SUBST($2_CFLAGS)

fi

# Set shell variable and define macro with test result
USE_$2=$GV_FOUND
if test $GV_FOUND = "yes" ; then
    AC_DEFINE([HAVE_$2],[],[Defined if we have $2])
fi

])

