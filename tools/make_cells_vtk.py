#!/usr/bin/env python3
"""A script to convert SWIFT cell dumps to XML for paraview visualization.

Example:

    $ python make_cells_vtk.py input_file output_file
"""
import argparse
import pandas as pd
import numpy as np


def make_cells_vtk(input_file, output_basename):
    """
    Convert cell dumps to XML file.

    This will ingest the input cell dump file and produce an XML file using
    the UnstructuredGrid prescription from VTK.

    Docs on the VTK format can be found here:
    https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html

    Args:
        input_file (str): The input file to convert.
        output_file (str): The output file to save the converted data.
    """

    # Open the cell dump csv with pandas (skipping the first non-header row
    # which only contains total particle counts)
    df = pd.read_csv(
        input_file,
        skiprows=[
            1,
        ],
    )

    # Get the cell data we want
    cell_locx = df.loc1
    cell_locy = df.loc2
    cell_locz = df.loc3
    cell_widthx = df.width1
    cell_widthy = df.width2
    cell_widthz = df.width3
    cell_gpart_count = df.gpart_count
    cell_hydro_count = df.hydro_count
    cell_stars_count = df.stars_count
    cell_mpi_rank = df.mpi_rank
    cell_id = df.name
    cell_super = df.super
    cell_hydro_super = df.hydro_super
    cell_hydro_hmax = df.hydro_h_max
    cell_stars_hmax = df.stars_h_max
    cell_depth = df.depth

    # Derive some extra flags
    cell_is_super = (cell_id == cell_super).astype(int)
    cell_is_hydro_super = (cell_id == cell_hydro_super).astype(int)

    # Count cells
    ncell = len(cell_locx)

    # Count points
    npoint = ncell * 8

    # Derive Points
    point_locx = np.zeros(npoint)
    point_locy = np.zeros(npoint)
    point_locz = np.zeros(npoint)
    for i in range(ncell):
        point_locx[i * 8 : i * 8 + 8] = np.array(
            [
                cell_locx[i] - cell_widthx[i] / 2,
                cell_locx[i] + cell_widthx[i] / 2,
                cell_locx[i] + cell_widthx[i] / 2,
                cell_locx[i] - cell_widthx[i] / 2,
                cell_locx[i] - cell_widthx[i] / 2,
                cell_locx[i] + cell_widthx[i] / 2,
                cell_locx[i] + cell_widthx[i] / 2,
                cell_locx[i] - cell_widthx[i] / 2,
            ]
        )
        point_locy[i * 8 : i * 8 + 8] = np.array(
            [
                cell_locy[i] - cell_widthy[i] / 2,
                cell_locy[i] - cell_widthy[i] / 2,
                cell_locy[i] + cell_widthy[i] / 2,
                cell_locy[i] + cell_widthy[i] / 2,
                cell_locy[i] - cell_widthy[i] / 2,
                cell_locy[i] - cell_widthy[i] / 2,
                cell_locy[i] + cell_widthy[i] / 2,
                cell_locy[i] + cell_widthy[i] / 2,
            ]
        )
        point_locz[i * 8 : i * 8 + 8] = np.array(
            [
                cell_locz[i] - cell_widthz[i] / 2,
                cell_locz[i] - cell_widthz[i] / 2,
                cell_locz[i] - cell_widthz[i] / 2,
                cell_locz[i] - cell_widthz[i] / 2,
                cell_locz[i] + cell_widthz[i] / 2,
                cell_locz[i] + cell_widthz[i] / 2,
                cell_locz[i] + cell_widthz[i] / 2,
                cell_locz[i] + cell_widthz[i] / 2,
            ]
        )

    # Derive PointData

    # Derive CellData

    # Write the VTK XML file
    with open(output_basename + ".vtu", "w") as f:
        f.write("<VTKFile type=UnstructuredGrid>\n")
        f.write("\t" + "<UnstructuredGrid>\n")
        f.write(
            "\t" * 2
            + f"<Piece NumberOfPoints={npoint} NumberOfCells={ncell}>\n"
        )
        f.write("\t" * 2 + "<Points>\n")
        f.write("\t" * 2 + "<DataArray type=Float32 NumberOfComponents=3>\n")
        f.write("\t" * 3 + " ".join(map(str, point_locx)) + "\n")
        f.write("\t" * 3 + " ".join(map(str, point_locy)) + "\n")
        f.write("\t" * 3 + " ".join(map(str, point_locz)) + "\n")
        f.write("\t" * 2 + "</DataArray>\n")
        f.write("\t" * 2 + "</Points>\n")
        f.write("\t" * 2 + "<Cells>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=connectivity>\n")
        f.write("\t" * 3 + " ".join(map(str, range(0, npoint))) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=offsets>\n")
        f.write("\t" * 3 + " ".join(map(str, range(0, npoint + 1, 8))) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 3 + "<DataArray type=UInt8 Name=types>\n")
        f.write("\t" * 3 + " ".join(map(str, [11] * ncell)) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 2 + "</Cells>\n")
        f.write("\t" * 2 + "<PointData>\n")
        f.write("\t" * 2 + "</PointData>\n")
        f.write("\t" * 2 + "<CellData>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=gpart_count>\n")
        f.write("\t" * 3 + " ".join(map(str, cell_gpart_count)) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=hydro_count>\n")
        f.write("\t" * 3 + " ".join(map(str, cell_hydro_count)) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=stars_count>\n")
        f.write("\t" * 3 + " ".join(map(str, cell_stars_count)) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=mpi_rank>\n")
        f.write("\t" * 3 + " ".join(map(str, cell_mpi_rank)) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=hydro_hmax>\n")
        f.write("\t" * 3 + " ".join(map(str, cell_hydro_hmax)) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=stars_hmax>\n")
        f.write("\t" * 3 + " ".join(map(str, cell_stars_hmax)) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=depth>\n")
        f.write("\t" * 3 + " ".join(map(str, cell_depth)) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=is_super>\n")
        f.write("\t" * 3 + " ".join(map(str, cell_is_super)) + "\n")
        f.write("\t" * 3 + "</DataArray>\n")
        f.write("\t" * 3 + "<DataArray type=Int32 Name=is_hydro_super>\n")
        f.write("\t" * 3 + " ".join(map(str, cell_is_hydro_super)) + "\n")
        f.write("\t" * 2 + "</DataArray>\n")
        f.write("\t" * 2 + "</CellData>\n")
        f.write("\t" * 2 + "</Piece>\n")
        f.write("\t" + "</UnstructuredGrid>\n")
        f.write("</VTKFile>\n")


if __name__ == "__main__":
    # Define command line arguments
    parser = argparse.ArgumentParser(
        description="Convert SWIFT cell dump files to XML "
        "files for paraview visualization."
    )
    parser.add_argument(
        "input_file", type=str, help="The input file to convert."
    )
    parser.add_argument(
        "output_basename",
        type=str,
        help="The output file to save the converted data (minus extension).",
    )
    args = parser.parse_args()

    # Run the main function
    make_cells_vtk(args.input_file, args.output_basename)
