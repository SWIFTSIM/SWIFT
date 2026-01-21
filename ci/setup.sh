#!/bin/bash

# This file is part of SWIFT.
# Copyright (C) 2026 p.w.draper@durham.ac.uk.
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

#+
#  Common setup and functions for jobs.
#
#  Peter W. Draper 20-JAN-2026.
#-

#  When exiting in error report current configuration options.
function ONEXIT {
   if test "$?" != 0; then
      echo "Current configuration: $(grep "\./configure" config.log)"
   fi
}
trap ONEXIT EXIT

#  Wrap calls to make so we can use -j <something> easily.
#  Also suppress output from libtool. Those logs are too big for gitlab.
function do_make {
    make -j 4 V=0 $*
}

#  Wrap ./configure calls so we can echo that line to the log.
function do_configure {
    echo "## CONFIGURE: $*"
    ./configure $*
}

#  Run a command and only show the output if the command fails, so if a test
#  fails, otherwise we only really need to see success.
#  Captures output to a temporary file.
function do_run {
    local tmp
    tmp=$(mktemp) || return 1

    if ! "$@" >"$tmp" 2>&1; then
        echo "Command failed: $*" >&2
        echo "Output:" >&2
        tail -2000 "$tmp" >&2
        rm -f "$tmp"
        exit 1
    fi

    rm -f "$tmp"
}

# More chat from the scripts. Re-enable if desperate.
#set -e
#set -x
