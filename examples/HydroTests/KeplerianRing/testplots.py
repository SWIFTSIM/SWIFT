"""
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2017
#
# Josh Borrow (joshua.borrow@durham.ac.uk)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# -----------------------------------------------------------------------------
#
# This program creates test plots for the initial condition generator provided
# for the Keplerian Ring example.
#
###############################################################################
"""


if __name__ == "__main__":
    import yt

    data = yt.load("initial_conditions.hdf5")

    plot = yt.ParticlePlot(
        data,
        "particle_position_x",
        "particle_position_y",
        "density"
    )

    plot.save("test.png")


