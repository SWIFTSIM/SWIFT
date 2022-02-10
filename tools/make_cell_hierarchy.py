#!/usr/bin/env python3
"""
Usage:
    make_cell_hierarchy.py --input <LIST OF FILES> --output <OUTPUT PREFIX>
                           [--serve]

Where <LIST OF FILES> is a list of input files (or file pattern matching that
list) corresponding to the cell hierarchy on different MPI ranks for a single
time step.
<OUTPUT PREFIX> is the name of the output .html and .csv files that will be
created. If the prefix path contains a directory that does not exist, this
directory will be created.
If the optional argument --serve is provided, the script will start a local
server that can run the generated .html page and opens it in a browser.

Based on a bash script written by Loic Hausammann.

This file is part of SWIFT.

Copyright (C) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

import argparse
import os
import sys

# Get the script folder, we need to retrieve the html template from it
scriptfolder = os.path.dirname(sys.argv[0])

# Handle the command line.
parser = argparse.ArgumentParser(description="Create cell hierarchy web page.")

parser.add_argument("input", help="Input file(s)", nargs="+")
parser.add_argument("output", help="Output file prefix")
parser.add_argument(
    "--serve",
    help="Run a local web server and display the generated web page?",
    action="store_true",
)
args = parser.parse_args()

# First check if the output prefix contains a folder and create it to make
# sure it exists
folder, fileprefix = os.path.split(args.output)
if len(folder) > 0:
    os.makedirs(folder, exist_ok=True)
else:
    folder = None

# Accumulate all input files into a single file
with open("{0}.csv".format(args.output), "w") as ofile:
    for fname in sorted(args.input):
        if not os.path.exists(fname):
            print('Error: "{0}" does not exist!'.format(fname))
            exit(1)
        with open(fname, "r") as ifile:
            ofile.write(ifile.read())

# Read the html template
with open("{0}/data/cell_hierarchy.html".format(scriptfolder), "r") as ifile:
    html = ifile.read()

# Replace the old csv file with the actual csv file
csvname = "{0}.csv".format(fileprefix)
html = html.replace("cell_hierarchy.csv", csvname)

# Write out the new html file
with open("{0}.html".format(args.output), "w") as ofile:
    ofile.write(html)

if args.serve:
    import http.server
    import socketserver
    import webbrowser

    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=folder, **kwargs)

    print("Running a local server. Use CTRL+C to stop it.")
    PORT = 8000
    found_port = False
    while not found_port:
        try:
            with socketserver.TCPServer(("", PORT), Handler) as httpd:
                found_port = True
                print("serving at port", PORT)
                webbrowser.open(
                    "http://localhost:{0}/{1}.html".format(PORT, fileprefix)
                )
                httpd.serve_forever()
        except:
            PORT += 1
