from argparse import ArgumentParser, RawTextHelpFormatter


class OptionsParser:
    def __init__(self):
        self.parser = ArgumentParser(
            description="Generate a cube of particles with a turbulent velocity field.",
            formatter_class=RawTextHelpFormatter,
        )

        self.parser.add_argument(
            "-Nside",
            "-nside",
            dest="nside",
            type=int,
            help="Number of particles per side of the box (required).",
            required=True,
        )

        self.parser.add_argument(
            "-o",
            metavar="outfile",
            dest="outfile",
            help="Name of output file.\n" + " [Default = ics_cloud.dat]",
            default="ics_cloud.dat",
        )

        self.parser.add_argument(
            "-format",
            dest="format",
            type=int,
            help="Format of output file:\n0 = ASCII,\n"
            + "1 = Gadget binary format 1,\n"
            + "2 = Gadget binary format 2,\n"
            + "3 = HDF5 format. \n"
            + " [Default = 3]",
            default=3,
        )

        self.parser.add_argument(
            "-m",
            "-mass",
            dest="mass",
            type=float,
            help="Total gas mass (in 10^10 Msun).\n" + " [Default = 1]",
            default=1.0,
        )

        self.parser.add_argument(
            "-npow",
            dest="npow",
            type=float,
            help="Power index of the power spectrum.\n" + " [Default = -4]",
            default=-4.0,
        )

        self.parser.add_argument(
            "-ngrid",
            dest="ngrid",
            type=int,
            help="Number of grid points per dimension for\n"
            + "the turbulent velocity field.\n"
            + " [Default = 256]",
            default=256,
        )

        self.parser.add_argument(
            "-l",
            "-lside",
            dest="lside",
            type=float,
            help="length side of the cube (in Mpc).\n" + " [Default = 1]",
            default=1.0,
        )

        self.parser.add_argument(
            "-a",
            "-alpha",
            dest="alpha",
            type=float,
            help="Ratio of turbulent energy to the\n"
            + "magnitude of gravitational energy.\n"
            + " [Default = 0.5]",
            default=0.5,
        )

    #        self.parser.add_argument("-b", "-beta",
    #                            dest     = "beta",
    #                            type     = float,
    #                            help     = "Ratio of rotational energy to the\n"+\
    #                                       "magnitude of gravitational energy.\n"+\
    #                                       " [Default = None]",
    #                            default  = -1)

    #        self.parser.add_argument("--units",
    #                            dest     = "units",
    #                            help     = "Change units to Msol/Parsec/km s^{-1} ",
    #                            action   = "store_true")

    def get_args(self):
        return self.parser.parse_args()
