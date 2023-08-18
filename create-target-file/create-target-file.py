import os
cwd = os.getcwd()

print("Current working directory:", cwd)


nodes_to_use = -1
mpi_tasks = -1
cpus_per_task = -1
tasks_per_core = -1
tasks_per_node = -1
cluster_architecture = "unknown"
processor_count = -1   # for -np=


# Sanity checks:
#  - if a [cwd] directive is used, create target.tgt in there;
#    otherwise, create it in $PWD, i..e in current directory
#  - if no cluster arch is specified, default to [cosma8]
#  - 
