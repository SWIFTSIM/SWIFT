import re
import os
import argparse
import math
import numpy as np
from scipy.interpolate import RegularGridInterpolator


def parse_interpolation_data(line, verbose):
    """
    Extract the parameters and the values of energy from the interpolation message.
    """
    pattern = re.compile(
        r"""
        interp->Nx=(?P<Nx>\d+)\s+
        interp->Ny=(?P<Ny>\d+)\s+
        interp->xmin=(?P<xmin>-?\d+\.?\d*)\s+
        interp->ymin=(?P<ymin>-?\d+\.?\d*)\s+
        interp->dx=(?P<dx>-?\d+\.?\d*)\s+
        interp->dy=(?P<dy>-?\d+\.?\d*)\s+
        idx=(?P<idx>-?\d+)\s+
        idy=(?P<idy>-?\d+)\s+
        out_of_boundary_type=(?P<boundary>none|zero|const)
        ( # If out of boundary const
            \s+cell_to_get=(?P<cell_to_get>\d+)
            \s+E\[(?P<row>\d+).*\]\[(?P<col>\d+).*\]=(?P<E>-?\d+\.?\d*[eE]?[+-]?\d*)
        )?
        ( # if inside boundaries
            \s+E\[idx\]\[idy\]=(?P<E00>-?\d+\.?\d*[eE]?[+-]?\d*)
            \s+E\[idx\]\[idy\+1\]=(?P<E01>-?\d+\.?\d*[eE]?[+-]?\d*)
            \s+E\[idx\+1\]\[idy\]=(?P<E10>-?\d+\.?\d*[eE]?[+-]?\d*)
            \s+E\[idx\+1\]\[idy\+1\]=(?P<E11>-?\d+\.?\d*[eE]?[+-]?\d*)
        )?
    """,
        re.VERBOSE,
    )

    match = pattern.search(line)

    if not match:
        return None

    data = match.groupdict()

    # Expected types for each data
    expected_types = {
        "Nx": int,
        "Ny": int,
        "xmin": float,
        "ymin": float,
        "dx": float,
        "dy": float,
        "idx": int,
        "idy": int,
        "E00": float,
        "E01": float,
        "E10": float,
        "E11": float,
        "cell_to_get": int,
        "row": int,
        "col": int,
        "E": float,
    }

    # Convert to appropriate type
    for key, value in data.items():
        if value is not None and key in expected_types:
            try:
                data[key] = expected_types[key](value)
            except (ValueError, TypeError):
                # To get specific case where it does not works
                if verbose:
                    print(
                        f"Warning: Could not convert '{key}' with value '{value}' to {expected_types[key].__name__}."
                    )
                pass

    if verbose:
        print(f"Interpolation's properties : {data}")

    return data


def parse_stellar_properties(line, verbose):
    """
    Extract the stellar mass, metallicity and the given preSN.
    """
    pattern = re.compile(
        r"Star_type=(?P<type>continuous|single).*"
        r"init_mass\[M_odot\]=(?P<mass>-?\d+\.?\d*[eE]?[+-]?\d*).*"
        r"metallicity\[Z_odot\]=(?P<metallicity>-?\d+\.?\d*[eE]?[+-]?\d*).*"
        r"Energy\[erg/yr\]=(?P<energy>-?\d+\.?\d*[eE]?[+-]?\d*)"
    )
    match = pattern.search(line)
    if match:
        # Extracing stellar data
        data = match.groupdict()
        data["mass"] = float(data["mass"])
        data["metallicity"] = float(data["metallicity"])
        data["energy"] = float(data["energy"])

        if verbose:
            print(f"Stellar properties : {data}")
        return data
    return None


def interpolate_2d(params, mass, metallicity, verbose, num_error):
    """
    Apply an interpolation from the given parameters and using scipy.
    """

    # Coordinate in the mass-metallicity log space
    # TODO: Do the case of negative values
    log_mass = math.log10(mass) if mass > 0 else float("-inf")
    log_metallicity = math.log10(metallicity) if metallicity > 0 else float("-inf")

    # Recalculate the grid values from the interpolation's properties
    x = (log_metallicity - params["xmin"]) / params["dx"]
    y = (log_mass - params["ymin"]) / params["dy"]
    idx = int(x)
    idy = int(y)
    fx = x - idx
    fy = y - idy

    # Check out of boundaries condition
    if x < params["Nx"] - 1 and y < params["Ny"] - 1 and x >= 0 and y >= 0:
        # calculate the 4 corner of the grid's targeted cell
        x0 = 10 ** (params["xmin"] + params["dx"] * idx)
        x1 = 10 ** (params["xmin"] + params["dx"] * (idx + 1))
        y0 = 10 ** (params["ymin"] + params["dy"] * idy)
        y1 = 10 ** (params["ymin"] + params["dy"] * (idy + 1))
        if verbose:
            print(
                f"x={x}, y={y}, idy={idx}, idy={idy}, fx={fx}, fy={fy}, x0={x0}, x1={x1}, y0={y0}, y1={y1}"
            )
    else:
        print("Out of boundaries")
        if params["boundary"] == "const":
            row = min(max(idx, 0), params["Nx"] - 1)
            col = min(max(idy, 0), params["Ny"] - 1)
            cell_to_get = row * params["Ny"] + col
            print(
                f"x={x}, y={y}, idy={idx}, idy={idy}, fx={fx}, fy={fy}, row={row}, col={col}, cell_to_get={cell_to_get}"
            )
            if (
                idx != params["idx"]
                or idy != params["idy"]
                or row != params["row"]
                or col != params["col"]
                or cell_to_get != params["cell_to_get"]
            ):
                num_error[0] += 1
                print(
                    f"Error : the indices found are in opposition with the SWIFT's ones"
                )
            return [None]

        elif params["boundary"] == "zero":
            print(
                f"x={x}, y={y}, idy={idx}, idy={idy}, fx={fx}, fy={fy}, row={row}, col={col}, cell_to_get={cell_to_get}"
            )
            if idx != params["idx"] or idy != params["idy"]:
                num_error[0] += 1
                print(
                    f"Error : the indices found are in opposition with the SWIFT's ones"
                )
            return [(0)]

        else:
            print("Error : type of boundary not found or not expected")
            return None

    # Interpolation by scipy
    x_coords = np.array([x0, x1])
    y_coords = np.array([y0, y1])
    energy_grid = np.array(
        [[params["E00"], params["E01"]], [params["E10"], params["E11"]]]
    )

    # By default bilinear
    f = RegularGridInterpolator((x_coords, y_coords), energy_grid)
    E_interp = f(np.array([metallicity, mass]))

    return E_interp


def main():
    parser = argparse.ArgumentParser(description="Test the 2D interpolation in swift")
    parser.add_argument(
        "log_file", type=str, nargs="+", help="path to the output.log file."
    )

    # Option for getting more details
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print more information about each individual file's calculation.",
    )

    args = parser.parse_args()

    num_error = [(0)]

    for file in args.log_file:

        if os.path.exists(file):
            print(f"Verification ocurring on file '{file}'")
        else:
            print(f"Error : The file '{file}' do not exist.\n")
            num_error += 1
            continue

        interp_data = None
        stellar_data = None

        with open(file, "r") as f:
            # Do we want to have the same line read many times ?
            # Knowing that the first time it appear correspond to the true calculation it should be ok to let it like this
            interp = False
            stellar = False
            for line in f:
                if "interpolate_2d: interp->Nx" in line and not interp:
                    interp_data = parse_interpolation_data(line, args.verbose)
                    interp = True
                elif (
                    "stellar_evolution_compute_preSN_properties: Star_type" in line
                    and not stellar
                ):
                    stellar_data = parse_stellar_properties(line, args.verbose)
                    stellar = True

                # The first time both lines are found is ok
                if interp_data and stellar_data:
                    break

        if not interp_data:
            print("Error : Interpolation's data are not found in this file.\n")
            num_error[0] += 1
            continue

        if not stellar_data:
            print("Error : Stellar data are not found in this file.\n")
            num_error += 1
            continue

        # We use the data we got to interpolate
        interpolated_energy = interpolate_2d(
            interp_data,
            stellar_data["mass"],
            stellar_data["metallicity"],
            args.verbose,
            num_error,
        )

        interpolated_energy = interpolated_energy[0]

        actual_energy = stellar_data["energy"]

        # Calculation of the error
        if interpolated_energy != None:
            if actual_energy != 0:
                error = (
                    abs(interpolated_energy - actual_energy) / abs(actual_energy) * 100
                )
            else:
                error = abs(interpolated_energy - actual_energy)

        # Print the results
        print("--- Interpolation's verification ---")
        print(f"The values of energy retrieved in Swift : {actual_energy:e} erg/yr")
        if interpolated_energy != None:
            print(f"The values interpolated : {interpolated_energy:e} erg/yr")
            print(
                f"Absolute error : {abs(interpolated_energy - actual_energy):e} erg/yr"
            )
            print(f"Relative error : {error:.4f}%")

            if error < 1:  # tolerance threshold
                print("The interpolation is correct.")
            else:
                print("The error on the interpolation is too big.")
                num_error[0] += 1

        print("\n ")

    print(f"In total the number of errors is {num_error[0]}")


if __name__ == "__main__":
    main()
