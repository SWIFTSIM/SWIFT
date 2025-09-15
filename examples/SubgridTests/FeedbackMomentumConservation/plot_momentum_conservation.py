################################################################################
# This file is part of SWIFT.
# Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_data(file_path):
    """
    Load data from a specified file path into a pandas DataFrame.

    This function reads a plain text data file, assuming it has no header,
    and converts it into a pandas DataFrame with predefined column names.

    Parameters
    ----------
    file_path : str
        The full or relative path to the data file (e.g., 'statistics.txt').

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the loaded data. Returns None if the file is
        not found or another error occurs during loading.

    Raises
    ------
    FileNotFoundError
        If the file specified by `file_path` does not exist.
    Exception
        For other potential errors during the file loading or processing.
    """
    try:
        # Define the column names for the dataset
        columns = [
            "Step",
            "Time",
            "Scale Factor",
            "Redshift",
            "Total Mass",
            "Gas Mass",
            "Dark Matter Mass",
            "Sink Mass",
            "Stellar Mass",
            "Black Hole Mass",
            "Gas Metal Mass",
            "Star Metal Mass",
            "BH Metal Mass",
            "Kinetic Energy",
            "Thermal Energy",
            "Potential Energy",
            "Radiated Energy",
            "Gas Entropy",
            "CoM_x",
            "CoM_y",
            "CoM_z",
            "Momentum_x",
            "Momentum_y",
            "Momentum_z",
            "AngMom_x",
            "AngMom_y",
            "AngMom_z",
            "BH Accretion Rate",
            "BH Accreted Mass",
            "BH Subgrid Mass",
            "Total H Mass",
            "Molecular H Mass",
            "Atomic H Mass",
            "Total He Mass",
            "Magnetic Energy",
            "DivB Error",
            "Cross Helicity",
            "Magnetic Helicity",
            "BH Luminosity",
            "BH Jet Power",
        ]

        # Load the data, ignoring any header comments
        data = np.loadtxt(file_path)

        # Convert the NumPy array to a pandas DataFrame
        df = pd.DataFrame(data, columns=columns)

        return df

    except FileNotFoundError:
        print(f"Error: The file at {file_path} was not found.")
        return None
    except Exception as e:
        print(f"An error occurred while loading the data: {e}")
        return None


# Path to the statistics.txt file
file_path = "./statistics.txt"

# Load the data
df = load_data(file_path)

# Calculate the total linear momentum
df["Momentum_Total"] = np.sqrt(
    df["Momentum_x"] ** 2 + df["Momentum_y"] ** 2 + df["Momentum_z"] ** 2
)
delta_p = np.abs(np.diff(df["Momentum_Total"]))

# Time to plot!
fig, axes = plt.subplots(figsize=(11, 4.8), nrows=1, ncols=2)

# Plot the total momentum
axes[0].plot(df["Time"], df["Momentum_Total"], label="Total Momentum", color="blue")
axes[0].set_xlabel("Time [Gyr]")
axes[0].set_ylabel(r"Momentum [$10^{10} M_\odot$ km/s]")
axes[0].set_yscale("log")
axes[0].grid(True)
axes[0].legend()

# Plot the change in total momentum (delta p) between two snapshots
axes[1].plot(df["Time"][1:], delta_p, label="Delta Total Momentum", color="red")
axes[1].set_xlabel("Time [Gyr]")
axes[1].set_yscale("log")
axes[1].grid(True)
axes[1].legend()

plt.show()
