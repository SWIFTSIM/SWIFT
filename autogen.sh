#! /bin/sh

#  Update generated configuration files, i.e. do work so that a
#  developer checkout can be configured.

./tools/update-modules
autoreconf --install --symlink 

exit

