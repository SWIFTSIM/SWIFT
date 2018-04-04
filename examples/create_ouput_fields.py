#!/usr/bin/env python3

from optparse import OptionParser
from regex import search
from os import path, remove

# fields to do
todo = [
    "chemistry", "cooling",
    "gravity", "stars",
    "hydro"]

# number of particle type
N_type = 6


def parseOptions():
    """
    Parse options.

    Returns
    -------

    opt:
        Option structure

    args:
        Arguments
    """
    usage = """usage: %prog [options] args1
    If options are not provided, this script will read the config.h file \
provided (args1) and use theses values.
    """

    parser = OptionParser(usage=usage)

    parser.add_option("--with-cooling",
                      dest="cooling",
                      help="Cooling Type",
                      default="",
                      type=str,
                      metavar="STR")

    parser.add_option("--with-chemistry",
                      dest="chemistry",
                      help="Chemistry Type",
                      default="",
                      type=str,
                      metavar="STR")

    parser.add_option("--with-hydro",
                      dest="hydro",
                      help="Hydro Type",
                      default="",
                      type=str,
                      metavar="STR")

    parser.add_option("--with-stars",
                      dest="stars",
                      help="Stars Type",
                      default="default",
                      type=str,
                      metavar="STR")

    parser.add_option("--with-gravity",
                      dest="gravity",
                      help="Gravity Type",
                      default="default",
                      type=str,
                      metavar="STR")

    parser.add_option("-o", "--output",
                      dest="o",
                      help="Output filename",
                      type=str,
                      default=None,
                      metavar="STR")

    parser.add_option("-d", "--example_directory",
                      dest="d",
                      help="Directory containing the examples",
                      type=str,
                      default=None,
                      metavar="STR")

    return parser.parse_args()


def readConfig(filename):
    """
    Read the config.h file and parse it to keep
    only the definitions.

    Parameters
    ----------

    filename: str
        Filename to parse

    Returns
    -------

    options: list
        List of options
    """
    expression = "#define (.*) 1"
    options = []
    with open(filename, "r") as f:
        for line in f:
            tmp = search(expression, line)
            if tmp is not None:
                options.append(tmp.group(1))

    return options


def getConfig(opt, swift_options, current):
    """
    Parse the config.h options and extract the required
    parameter. If not able to find an option, return False

    Parameters
    ----------

    opt:
        Option structure

    swift_options: list
        Option in config.h

    current: str
        Current option

    Returns
    -------

    option: str
        name of the option (in lower case)
    """
    # if given in options, returns
    if getattr(opt, current) != "":
        return getattr(opt, current)

    # get current option
    caps = current.upper()
    expression = "%s_(.*)" % caps
    if current == "hydro":
        expression = "(.*)_SPH"
    for i in swift_options:
        tmp = search(expression, i)
        if tmp is not None:
            return tmp.group(1).lower()
    return False


def readYamlFile(filename):
    """
    Read quickly a yaml file.

    Parameters
    ----------

    filename: str
        file to read

    Returns
    -------

    d: dict
        key -> section name, value -> list of parameters
    """
    d = {}
    key = None
    section = "(PartType\d+):"
    with open(filename, "r") as f:
        for line in f:
            # check for section
            tmp = search(section, line)
            if tmp is not None:
                key = tmp.group(1)
                d[key] = []
                continue

            # skip if not a parameter or outside
            # a section
            if ":" not in line or key is None:
                continue

            # ensure end of line
            if "\n" not in line:
                line += "\n"

            # write line
            d[key].append(line)

    return d


def generateFile(opt):
    """
    Generate an output fields file from examples.

    Parameters
    ----------

    opt:
        Option structure

    Returns
    -------

    d: dict
        key -> section name, value -> list of parameters
    """
    d = {}
    # read all files
    for current in todo:
        filename = current + "_" + getattr(opt, current) + ".yml"
        filename = opt.d + "/" + filename
        data = readYamlFile(filename)

        for tpe in range(N_type):
            name = "PartType%i" % tpe
            if name not in data:
                continue

            if name not in d:
                d[name] = data[name]

            else:
                d[name].extend(data[name])
    return d


def writeOutputFields(d, opt):
    """
    Write the output.

    Parameters
    ----------

    d: dict
        key -> section name, value -> list of parameters

    opt:
        Option structure

    """
    with open(opt.o, "w") as f:
        for tpe in range(N_type):
            name = "PartType%i" % tpe
            if name not in d:
                continue
            f.write("%s:\n" % name)
            for i in d[name]:
                f.write(i)

            f.write("\n")


if __name__ == "__main__":
    # parse option
    opt, args = parseOptions()

    # check inputs
    # arguments
    if len(args) != 1:
        raise IOError("A file should be provided")
    filename = args[0]

    # output file
    if opt.o is None:
        raise IOError("The output option '-o' is required")

    # example directory
    if opt.d is None:
        raise IOError("The example directory option '-d' is required")

    # read configuration file
    swift_options = readConfig(filename)

    # get correct configuration
    for current in todo:
        tmp = getConfig(opt, swift_options, current)
        if tmp is False:
            raise IOError("Unable to get field %s" % current)
        setattr(opt, current, tmp)

    # generate and write output_fields file
    d = generateFile(opt)
    writeOutputFields(d, opt)
