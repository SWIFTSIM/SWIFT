import os
import shutil
import glob
import swiftsimio

the_folder = "./"
snapshot_key = "Monopole_"
parameter_file = "Monopole.yml"

# getting snapshot addresses
addressbook = glob.glob(the_folder + "*.hdf5")
snapshots = sorted(
    [
        addressbook[i]
        for i in range(len(addressbook))
        if snapshot_key in os.path.splitext(addressbook[i])[0]
    ]
)

# executing plotSolution.py for all snapshots
for i in range(len(snapshots)):
    print("")
    print(f"Processing snapshot {i}")
    filename = os.path.splitext(snapshots[i])[0]
    ext = ".png"
    file_to_open = filename + ext
    cwd = os.path.join(os.getcwd(), "div_maps.py " + snapshots[i] + " " + file_to_open)
    os.system("{} {}".format("python3", cwd))
