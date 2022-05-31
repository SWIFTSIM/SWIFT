#!/usr/bin/env python3
"""
Usage:
    recover_restart_parameters.py <STEP> <LIST OF FILES>

where <STEP> is a simulation time step number, and <LIST OF FILES> are all
the 'used_parameters.yml[.stepno]' files that were produced by the run.

This script will reconstruct the parameter file as it was used by SWIFT when
it (last) ran time step <STEP>, and will take into account any potential
parameter changes that happened before restarting.

Copyright (C) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
All Rights Reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import yaml
import argparse
import numpy as np
import datetime
import re
import sys

argparser = argparse.ArgumentParser(
    "Recover the parameters used for a specific time step."
)
argparser.add_argument(
    "step", type=int, help="Step for which we want to recover the parameters."
)
argparser.add_argument("file", nargs="+", help="used_parameters.yml* file(s) to parse.")
args = argparser.parse_args()


def get_step(filename):
    """
    Extract the step number from a used_parameters.yml.stepno file.
    Return 0 if no step number was attached (corresponding to the first step).
    """
    match = re.search("used_parameters.yml.(\d+)", filename)
    if match is None:
        return 0
    else:
        return int(match[1])


# make a data array that contains the timestamp, step number and filename
# of each input file
files = np.zeros(
    len(args.file),
    dtype=[("timestamp", np.uint64), ("step", np.uint32), ("filename", "U100")],
)
for i, file in enumerate(args.file):
    # make sure we are parsing a used_parameters file
    if not "used_parameters.yml" in file:
        raise ValueError(f'Incompatible filename: "{file}"!')
    # store the filename and step number
    files[i]["filename"] = file
    files[i]["step"] = get_step(file)
    # extract the time stamp from the file header
    with open(file, "r") as ifile:
        for line in ifile.readlines():
            if "current date:" in line:
                timestr = " ".join(line.split()[3:])
                date = datetime.datetime.strptime(timestr, "%H:%M:%S %Y-%m-%d %Z")
                files[i]["timestamp"] = date.timestamp()

# sort the files according to timestamp
isort = np.argsort(files["timestamp"])
files = files[isort]

# get the index of the first step
ifirst = np.argmax(files["step"] == 0)
# discard any step files that are older
files = files[ifirst:]

# filter out all the steps that contributed to the step we want
mask = files["step"] <= args.step
files = files[mask]
# only the first (step 0) and last file matter
files = files[[0, -1]]


def update_dictionary(d, nd):
    """
    Recursively update the contents of dictionary 'd' with that of dictionary
    'nd'.
    """
    for key in nd:
        if key in d and isinstance(nd[key], dict):
            update_dictionary(d[key], nd[key])
        else:
            d[key] = nd[key]


# now create the final parameter file
params = {}
for file in files["filename"]:
    with open(file, "r") as ifile:
        this_params = yaml.safe_load(ifile)
    # if no parameters changed compared to step 0, this_params is empty
    if not this_params is None:
        update_dictionary(params, this_params)

# dump the result to the stdout
# note that pyYAML will format array parameters as lists instead of inline
# arrays. There is not much we can do about that.
yaml.safe_dump(params, stream=sys.stdout, default_flow_style=False)
