#!/bin/bash
#
# Usage:
#  process_memuse_logs_MPI rank nprocess
#
# Description:
#  Process all the memuse report files in the current directory that
#  are output from the given rank.
#  Creating an analysis for each step and one for all the steps.
#
#  The input files are created by a run configured for memuse reporting
#  (--enable-memuse-reports) should be named "memuse_report-rank<n>-step<m>.dat"
#  in the current directory.
#
#  All located files will be processed using "nprocess" concurrent
#  processes. The output for each step will be named memuse_report-rank<n>-step<m>.log
#  and the overall analysis will be called memuse_report-all-rank<n>.log.
#
# This file is part of SWIFT:
#
#  Copyright (C) 2019 Peter W. Draper (p.w.draper@durham.ac.uk)
#  All Rights Reserved.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#  Handle command-line
if test "$2" == ""; then
    echo "Usage: $0 rank nprocess"
    exit 1
fi
RANK=$1
NPROCS=$2

#  Locate script.
SCRIPTHOME=$(dirname "$0")

#  Find all report files. Use version sort to get into correct order.
files=$(ls -v memuse_report-rank${RANK}-step*.dat)
if test $? != 0; then
    echo "Failed to find any memuse report files"
    exit 1
fi

#  Construct list of input and output names.
list=""
for f in $files; do
    output=$(echo $f| sed 's,.dat,.log,')
    list="$list $f $output"
done

#  And process them.
echo "Processing memuse report files..."
echo $list | xargs -P $NPROCS -n 2 /bin/bash -c "${SCRIPTHOME}/analyse_memuse_logs.py \$0 > \$1"

#  Now process the overall file, if more than one file given. 
n=$(echo $list| wc -w)
if test $n -gt 2; then
    echo "Processing _all_ memuse report files..."
    ${SCRIPTHOME}/analyse_memuse_logs.py $files > memuse_report-all-rank${RANK}.log
fi

echo "Finished"

exit
