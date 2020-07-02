#!/bin/bash

#-------------------------------------------------------------
# Don't forget to adapt the program_suffix appropriately!
#-------------------------------------------------------------


DEFAULTFLAGS=""
DEFAULTFLAGS="$DEFAULTFLAGS"" --enable-mpi=no"
DEFAULTFLAGS="$DEFAULTFLAGS"" --disable-doxygen-doc"


DEBUGFLAGS=""
# DEBUGFLAGS="$DEBUGFLAGS"" --enable-debug"
# DEBUGFLAGS="$DEBUGFLAGS"" --enable-sanitizer"
# DEBUGFLAGS="$DEBUGFLAGS"" --enable-undefined-sanitizer"
# DEBUGFLAGS="$DEBUGFLAGS"" --enable-debugging-checks"
# DEBUGFLAGS="$DEBUGFLAGS"" --enable-optimization=no"
# DEBUGFLAGS="$DEBUGFLAGS"" --enable-task-debugging"


# without Ivanova
GIZMOFLAGS=""
GIZMOFLAGS="$GIZMOFLAGS""--with-hydro-dimension=2"
GIZMOFLAGS="$GIZMOFLAGS"" --with-hydro=gizmo-mfv"
# GIZMOFLAGS="$GIZMOFLAGS"" --with-riemann-solver=hllc"
GIZMOFLAGS="$GIZMOFLAGS"" --with-riemann-solver=exact"
GIZMOFLAGS="$GIZMOFLAGS"" --enable-ivanova-surfaces"


LIBFLAGS=""
LIBFLAGS="$LIBFLAGS"" --with-parmetis"
LIBFLAGS="$LIBFLAGS"" --with-jemalloc"
LIBFLAGS="$LIBFLAGS"" --with-hdf5=$HDF5_ROOT/bin/h5pcc"

EXTRA_CFLAGS=""
# EXTRA_CFLAGS="$EXTRA_CFLAGS"" -save-temps"


ADDFLAGS=""
# ADDFLAGS="$ADDFLAGS"" --with-kernel=wendland-C6"
# ADDFLAGS="$ADDFLAGS"" --disable-vec"
# ADDFLAGS="$ADDFLAGS"" --disable-hand-vec"



allflags="$LIBFLAGS ""$GIZMOFLAGS ""$DEFAULTFLAGS"" $ADDFLAGS"" $DEBUGFLAGS"" CFLAGS=$EXTRA_CFLAGS" 




if [ ! -f ./configure ]; then
    ./autogen.sh
fi








#---------------------------------
# generate exec filename suffix
#---------------------------------

program_suffix=-"2d"
# program_suffix=-"2d-debug"
# program_suffix=-"2d-debug-testing" # use "testing" additional suffix when you manually tinkered with swift's insides
# program_suffix=-"2d-clean"
# program_suffix=-"2d-clean-debug"
# program_suffix=-"2d-testing"





# make clean
# ./configure $allflags
make -j | tee -a compile_log





#--------------------------------------
# store what this compilation was
#--------------------------------------
echo "$allflags" | tr -d \\n | sed -r 's/\s+/ /g' > .last_compile
echo >> .last_compile # tr -d \\n removes all newlines, including the last one, so add one here


# rename executables
execname="./examples/swift""$program_suffix"
execname_mpi="./examples/swift""$program_suffix"-mpi
mv ./examples/swift "$execname"
echo "renamed ./examples/swift -> $execname"
if [ -f ./examples/swift_mpi ]; then 
    mv ./examples/swift_mpi "$execname_mpi"
    myecho "renamed ./examples/swift_mpi -> $execname_mpi"
fi
