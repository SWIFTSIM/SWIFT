#!/usr/bin/env python3
"""
This file generates a graphviz file that represents the SWIFT tasks
 dependencies.

Example: ./plot_task_dependencies.py dependency_graph_*.csv
"""
import sys
from pandas import read_csv
import numpy as np
from subprocess import call


def getGitVersion(f, git):
    """
    Read the git version from the file

    Parameters
    ----------

    f: str
        Filename

    git: str
        Git version of previous file

    Returns
    -------

    new_git: str
        Git version of current file
    """
    # read comment in csv file
    with open(f, "r") as f:
        line = f.readline()

    # check if really a comment
    if line[0] != "#":
        return None

    # remove trailing characters
    new_git = line[2:].rstrip()

    # check if previous and current are the same
    if git is not None and git != new_git:
        raise Exception("Files were not produced by the same version")

    return new_git


def appendSingleData(data0, datai):
    """
    Append two DataFrame together

    Parameters
    ----------

    data0: DataFrame
        One of the dataframe

    datai: DataFrame
        The second dataframe

    Returns
    -------

    data0: DataFrame
        The updated dataframe
    """

    # loop over all rows in datai
    for i, row in datai.iterrows():
        # get data
        ta = datai["task_in"][i]
        tb = datai["task_out"][i]
        ind = np.logical_and(data0["task_in"] == ta,
                             data0["task_out"] == tb)

        # check number of ta->tb
        N = np.sum(ind)
        if N > 1:
            raise Exception("Same dependency written multiple times %s->%s" %
                            (ta, tb))
        # if not present in data0
        if N == 0:
            data0.append(row)
        else:
            # otherwise just update the number of link
            ind = ind[ind].index[0]
            tmp = data0["number_link"][ind] + datai["number_link"][i]
            data0.at[ind, "number_link"] = tmp

    return data0


def appendData(data):
    """
    Append all the dataframe together

    Parameters
    ----------

    data: list
        List containing all the dataframe to append together

    Returns
    -------

    data: DataFrame
        The complete dataframe
    """
    N = len(data)
    if N == 1:
        return data[0]

    # add number link to data[0]
    for i in range(N-1):
        i += 1
        data[0] = appendSingleData(data[0], data[i])

    return data[0]


def taskIsStars(name):
    """
    Does the task concern stars?

    Parameters
    ----------

    name: str
        Task name
    """
    if "stars" in name or "spart" in name:
        return True
    return False


def taskIsHydro(name):
    """
    Does the task concern the hydro?

    Parameters
    ----------

    name: str
        Task name
    """
    if "_part" in name:
        return True
    if "density" in name and "stars" not in name:
        return True
    if "rho" in name:
        return True
    if "force" in name:
        return True
    if "xv" in name:
        return True

    task_name = [
        "sort",
        "ghost_in",
        "ghost",
        "ghost_out",
    ]
    if name in task_name:
        return True
    return False


def taskIsGravity(name):
    """
    Does the task concern the gravity?

    Parameters
    ----------

    name: str
        Task name
    """
    if "gpart" in name:
        return True
    if "grav" in name:
        return True
    return False


def writeTask(f, name, implicit, mpi):
    """
    Write the special task (e.g. implicit and mpi)

    Parameters
    ----------

    f: File
        File where to write the data

    name: str
        Task name

    implicit: int
        Is the task implicit

    mpi: int
        Is the task MPI related
    """
    # generate text
    txt = "\t " + name + "["

    if implicit:
        txt += "style=filled,fillcolor=lightgrey,"
    if mpi:
        txt += "shape=diamond,"

    if taskIsStars(name):
        txt += "color=darkorange1,"

    if taskIsHydro(name):
        txt += "color=blue3,"

    if taskIsGravity(name):
        txt += "color=red3,"

    # remove extra ','
    if txt[-1] == ",":
        txt = txt[:-1]
    txt += "];\n"

    # write it
    f.write(txt)


def writeHeader(f, data, git):
    """
    Write the header and the special tasks

    Parameters
    ----------

    f: File
        File where to write the data

    data: DataFrame
        The dataframe to write

    git: str
        The git version
    """
    # write header
    f.write("digraph task_dep {\n")
    f.write("\t # Header\n")
    f.write('\t label="Task dependencies for SWIFT %s";\n' % git)
    f.write("\t compound=true;\n")
    f.write("\t ratio=0.66;\n")
    f.write("\t node[nodesep=0.15];\n")

    f.write("\n")

    # write the special task
    f.write("\t # Special tasks\n")
    N = len(data)
    written = []
    # do task in
    for i in range(N):
        ta = data["task_in"][i]
        if ta in written:
            continue

        written.append(ta)
        writeTask(f, ta, data["implicit_in"][i], data["mpi_in"][i])

    # do task out
    for i in range(N):
        tb = data["task_out"][i]
        if tb in written:
            continue

        written.append(tb)
        writeTask(f, tb, data["implicit_out"][i], data["mpi_out"][i])

    f.write("\n")


def writeCluster(f, tasks, cluster):
    """
    Write a single cluster

    Parameters
    ----------

    f: File
        File where to write the data

    tasks: list
        List of all tasks in the cluster

    cluster: str
        Cluster name
    """
    f.write("\t subgraph cluster%s {\n" % cluster)
    f.write('\t\t label="";\n')
    for t in tasks:
        f.write("\t\t %s;\n" % t)
    f.write("\t };\n\n")


def writeClusters(f, data):
    """
    Write all the clusters

    Parameters
    ----------

    f: File
        File where to write the data

    data: DataFrame
        The dataframe to write
    """
    f.write("\t # Clusters\n")
    # get list of all the clusters
    clusters = data[["cluster_in", "cluster_out"]]
    clusters = np.unique(clusters)

    cluster_in = data["cluster_in"]
    cluster_out = data["cluster_out"]
    # loop over all clusters
    for cluster in clusters:
        # is it a cluster?
        if cluster == "None":
            continue

        # get all the task in current cluster
        ta = data["task_in"][cluster_in == cluster]
        tb = data["task_out"][cluster_out == cluster]

        # make them unique
        tasks = np.append(ta, tb)
        tasks = np.unique(tasks)

        # write current cluster
        writeCluster(f, tasks, cluster)

    f.write("\n")


def writeDependencies(f, data):
    """
    Write all the dependencies between tasks

    Parameters
    ----------

    f: File
        File where to write the data

    data: DataFrame
        The dataframe to write

    """
    f.write("\t # Dependencies\n")
    N = len(data)
    written = []
    max_rank = data["number_rank"].max()
    for i in range(N):
        # get data
        ta = data["task_in"][i]
        tb = data["task_out"][i]
        number_link = data["number_link"][i]

        # check if already done
        name = "%s_%s" % (ta, tb)
        if name in written:
            raise Exception("Found two same task dependencies")

        written.append(name)

        # write relation
        arrow = ""
        if data["number_rank"][i] != max_rank:
            arrow = ",style=dashed"
        f.write("\t %s->%s[label=%i%s]\n" %
                (ta, tb, number_link, arrow))


def writeFooter(f):
    """
    Write the footer

    Parameters
    ----------

    f: File
        File where to write the data
    """
    f.write("}")


if __name__ == "__main__":
    # get input
    filenames = sys.argv[1:]
    if len(filenames) < 1:
        raise Exception("You should provide at least a file name")

    # output
    dot_output = "dependency_graph.dot"
    png_output = "dependency_graph.png"

    # read files
    data = []
    git = None
    for f in filenames:
        tmp = read_csv(f, delimiter=",", comment="#")
        git = getGitVersion(f, git)
        data.append(tmp)

    data = appendData(data)

    # write output
    with open(dot_output, "w") as f:
        writeHeader(f, data, git)

        writeClusters(f, data)

        writeDependencies(f, data)

        writeFooter(f)

    call(["dot", "-Tpng", dot_output, "-o", png_output])

    print("You will find the graph in %s" % png_output)
