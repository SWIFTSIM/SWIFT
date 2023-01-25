import os

# Yes... this is overkill... but the main sphinx conf script looks for python
# scripts to execute... So we piggy back on this.
os.system("dot -Tpng full_dependency_graph_RT.dot -o full_dependency_graph_RT.png")
