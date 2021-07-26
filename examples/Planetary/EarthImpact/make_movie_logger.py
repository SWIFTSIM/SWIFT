#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import swiftsimio.visualisation.projection as vis
from copy import deepcopy
from shutil import copyfile
import os
from subprocess import call
sys.path.append("../../../csds/.libs/")
import libcsds as csds

boxsize = 80.
large_width = 15
small_width = 0.4
alpha_particle = 0

width = 0
basename = "index"
resolution = 1080
resolution_mini = resolution * 100 // 512
id_foc = 99839
v_max = 8.

traj = None


def selectParticles(pos, width):
    ind1 = np.logical_and(pos[:, 0] > 0, pos[:, 0] < width)
    ind2 = np.logical_and(pos[:, 1] > 0, pos[:, 1] < width)
    ind3 = np.logical_and(ind1, ind2)

    # avoid the zoom
    size = resolution_mini * width / resolution
    box = np.logical_and(pos[:, 0] > (width - size), pos[:, 1] < size)
    ind3 = np.logical_and(~box, ind3)

    return ind3


def getImage(parts, width, center, res):
    pos = parts["positions"]
    pos = pos - center + 0.5 * width
    pos /= width
    h = parts["smoothing_lengths"] / width

    m = np.ones(pos.shape[0])
    # Do the projection
    img = vis.scatter(pos[:, 0], pos[:, 1], m, h, res).T
    img /= width**3
    ind = img == 0
    img[ind] = img[~ind].min()
    img = np.log10(img)

    return img


def doProjection(parts, t, skip):
    plt.figure(figsize=(8, 8))
    global traj, id_foc

    # evolve in time the particles
    interp = csds.moveForwardInTime(basename, parts, t)

    # Check if some particles where removed
    ind = parts["smoothing_lengths"] == 0
    ind = np.arange(parts.shape[0])[ind]
    if len(ind) != 0:
        parts = np.delete(parts, ind)
        interp = np.delete(interp, ind)
        id_foc -= np.sum(ind < id_foc)

    # Get the position and the image center
    pos = interp["positions"]
    position = pos[id_foc, :]

    # save the trajectory
    if traj is None:
        traj = position[np.newaxis, :]
    else:
        traj = np.append(traj, position[np.newaxis, :], axis=0)

    # Do we want to generate the image?
    if skip:
        return parts

    # Compute the images
    img = getImage(interp, width, position, resolution)
    img[:resolution_mini, -resolution_mini:] = getImage(
        interp, 0.2 * boxsize, position, resolution_mini)
    box = width * np.array([0, 1, 0, 1])
    plt.imshow(img, origin="lower", extent=box, cmap="plasma",
               vmin=0, vmax=v_max)

    # plot the particles
    pos = interp["positions"] - position + 0.5 * width
    ind = selectParticles(pos, width)
    ms = 2
    plt.plot(pos[ind, 0], pos[ind, 1], "o", markersize=ms,
             alpha=alpha_particle, color="silver")
    plt.plot(pos[id_foc, 0], pos[id_foc, 1], "or", markersize=2*ms,
             alpha=alpha_particle)

    # plot time
    plt.text(0.5 * width, 0.95 * width,
             "Time = %0.2f Hours" % (t / (60 * 60)), color="w",
             horizontalalignment="center")

    # plot trajectory
    tr = deepcopy(traj) - position + 0.5 * width
    plt.plot(tr[:, 0], tr[:, 1], "-", color="w", alpha=alpha_particle)

    # plot scale
    plt.plot([0.05 * width, 0.15 * width], [0.05 * width, 0.05 * width],
             "-w")
    plt.text(0.1 * width, 0.06 * width, "%.2f R$_\oplus$" % (0.1 * width),
             horizontalalignment='center', color="w")

    # style
    plt.axis("off")
    plt.xlim(box[:2])
    plt.ylim(box[2:])
    plt.tight_layout()
    plt.style.use('dark_background')
    return parts


def skipImage(i):
    if os.path.isfile("output/image_%04i.png" % i):
        print("Skipping image %i" % i)
        return True
    else:
        return False


def doMovie(parts, t0, t1, N, init):
    times = np.linspace(t0, t1, N)
    for i, t in enumerate(times):
        print("Image %i / %i" % (i+1, N))
        skip = skipImage(i + init)
        parts = doProjection(parts, t, skip)
        if not skip:
            plt.savefig("output/image_%04i.png" % (i + init))
        plt.close()
    return init + N


def doZoom(parts, w_init, w_end, N, init, t0, t1, increase_alpha):
    global width, alpha_particle
    widths = np.linspace(w_init, w_end, N)
    alpha = np.linspace(-8, 0.2, N)
    if not increase_alpha:
        alpha = alpha[::-1]
    k = 5  # parameter for the steepness
    alpha = 1. / (1 + np.exp(- 2 * k * alpha))
    times = np.linspace(t0, t1, N)
    for i, w in enumerate(widths):
        print("Image %i / %i" % (i+1, N))
        skip = skipImage(i + init)
        width = w
        alpha_particle = alpha[i]
        parts = doProjection(parts, times[i], skip)
        if not skip:
            plt.savefig("output/image_%04i.png" % (i + init))
        plt.close()
    return init + N


def doStatic(init, number, ref):
    # copy the first picture
    for i in range(number):
        copyfile("output/image_%04i.png" % ref,
                 "output/image_%04i.png" % (i + init))
    return init + number


def doTitle(frames):
    plt.figure()
    plt.axis("off")
    box = np.array([0, 1, 0, 1])
    plt.xlim(box[:2])
    plt.ylim(box[2:])
    plt.tight_layout()
    plt.style.use('dark_background')

    style = {
        "verticalalignment": "center",
        "horizontalalignment": "center",
        "fontweight": "bold"
    }
    plt.text(0.5, 0.6, "Planetary Impact with the CSDS",
             **style)
    plt.text(0.5, 0.5, "L. Hausammann, J. Kegerreis and P. Gonnet 2020",
             **style)

    for i in range(frames):
        plt.savefig("output/image_%04i.png" % i)

    plt.close()
    return frames


def main():
    t0, t1 = csds.getTimeLimits(basename)
    parts = csds.loadSnapshotAtTime(basename, t0)

    # Do a few title frames
    init = doTitle(40)

    after_title = init
    frames_after_title = 10
    init += frames_after_title

    # do a zoom while moving forward in time
    init = doZoom(parts, large_width, small_width, 50,
                  init=init, t0=t0, t1=0.07 * t1,
                  increase_alpha=True)

    # Copy a few frames
    doStatic(after_title, frames_after_title,
             after_title + frames_after_title)
    init = doStatic(init, 10, init-1)
    # Do the main part of the movie
    init = doMovie(parts, 0.07 * t1, 0.15 * t1, N=250,
                   init=init)
    # copy a few frames
    init = doStatic(init, 10, init-1)

    # zoom out and finish the movie
    init = doZoom(parts, small_width, large_width, 607,
                  init=init, t0=0.15 * t1, t1=t1,
                  increase_alpha=False)

    # copy a few frames
    init = doStatic(init, 10, init-1)


if __name__ == "__main__":
    main()

    convert = "ffmpeg -i output/image_%04d.png -y -vcodec libx264 "
    convert += "-profile:v high444 -refs 16 -crf 0 "
    convert += "-preset ultrafast movie.mp4"
    call(convert, shell=True)
