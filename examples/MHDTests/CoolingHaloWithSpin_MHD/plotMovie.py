import subprocess
from glob import glob
import os


def import_all_snapshot_addr(folderaddr):
    addr_book = glob(folderaddr + "**/CoolingHalo_*.hdf5", recursive=True)
    addr_book = sorted(addr_book)
    image_names = [
        "image_" + addr.split(".hdf5")[0].split("_")[-1] + ".png" for addr in addr_book
    ]
    return addr_book, image_names


def plot_all_snapshots(folderaddr):
    snap_addr, output_addr = import_all_snapshot_addr(folderaddr)
    for saddr, iaddr in zip(snap_addr, output_addr):
        command = "python3 plotCorrelation.py " + saddr + " " + iaddr
        try:
            subprocess.call(command, shell=True)
        except subprocess.CalledProcessError as e:
            print(
                f"Encountered error number {e.returncode}, command output: {e.output}"
            )


def make_a_movie(framerate=1):
    command = f"ffmpeg -framerate {framerate} -i 'image_%04d.png' -avoid_negative_ts make_zero -r {framerate} -pix_fmt yuv420p video.mp4"
    try:
        subprocess.call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Encountered error number {e.returncode}, command output: {e.output}")


plot_all_snapshots("./CH_gb099172f_hyp=0.1_par=1.0_divv=0.0_smallHalo_1e4_correcteps")
# make_a_movie()
