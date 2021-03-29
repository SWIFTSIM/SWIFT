#! /bin/sh

#  Update generated configuration files, i.e. do work so that a
#  developer checkout can be configured.

if [ -f csds/Makefile.am ]; then
    fake_sub=0
else
    echo "Creating temporary (fake) submodule files"
    fake_sub=1
    mkdir csds/src csds/tests
    touch csds/Makefile.am csds/src/Makefile.am csds/tests/Makefile.am
fi

autoreconf --install --symlink

if [ $fake_sub -eq 1 ]; then
    echo "Removing fake submodule files"
    rm -rf csds/src csds/tests
    rm csds/Makefile.am csds/Makefile.in
fi

exit

