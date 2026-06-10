from __future__ import print_function

from sys import exc_info, exit
from os import path, remove
from numpy import array, zeros, uint32, int32, float32
from logging import warning
from struct import pack
from h5py import File

from libs.const import msol, parsec


def save_particles(ids, pos, vel, mass, u, outfile, format, boxsize, h):

    ngas = len(mass)

    if path.isfile(outfile):
        print("WARNING: File {} already exist.".format(outfile))
        print("Do yo want to overwrite it? Y/[N]")
        q = input()
        if q.lower() in ["y", "yes", "s", "si"]:
            remove(outfile)
        else:
            print("Exiting.")
            exit()

    if format == 0:
        # Openning file
        try:
            ofile = open(outfile, "w")
        except IOError as e:
            msg = "IO Error({0}): {1}".format(e.errno, e.strerror)
            warning(msg)
        except:
            print("Unexpected error: {}".format(exc_info()[0]))
            raise

        id_space = len("{}".format(ngas))

        # Preparing every line to print to the file
        for i in range(ngas):
            # Formatting particle attributes
            ie = "% d" % ids[i]
            me = "% 3.8e" % mass[i]
            rx = "% 3.8e" % pos[i][0]
            ry = "% 3.8e" % pos[i][1]
            rz = "% 3.8e" % pos[i][2]
            vx = "% 3.8e" % vel[i][0]
            vy = "% 3.8e" % vel[i][1]
            vz = "% 3.8e" % vel[i][2]
            ue = "% 3.8e" % u[i]

            # Right-align the strings
            outstring = "{0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(
                ie.rjust(id_space),
                me.rjust(12),
                rx.rjust(12),
                ry.rjust(12),
                rz.rjust(12),
                vx.rjust(12),
                vy.rjust(12),
                vz.rjust(12),
                ue.rjust(12),
            )
            # Write to file
            ofile.write(outstring)

        # Closing the file
        ofile.close()

    elif format == 1:
        npart = array([ngas, 0, 0, 0, 0, 0])
        Nmass = array([0, 0, 0, 0, 0, 0])
        # Linearizing the 3D-array of the position and velocity
        pos = pos.ravel()
        vel = vel.ravel()

        dummy = zeros(npart[0])
        time = 0.0  # double
        redshift = 0.0  # double
        flag_sfr = 0  # int
        flag_feedback = 0  # int
        bytesleft = 256 - 6 * 4 - 6 * 8 - 8 - 8 - 2 * 4 - 6 * 4
        fill = zeros(int(bytesleft / 4.0), dtype=int)  # int

        with open(outfile, "wb") as f:
            nbytes = 256
            # Header
            f.write(pack("i", nbytes))
            f.write(pack("i" * len(npart), *npart))
            f.write(pack("d" * len(Nmass), *Nmass))
            f.write(pack("d", time))
            f.write(pack("d", redshift))
            f.write(pack("i", flag_sfr))
            f.write(pack("i", flag_feedback))
            f.write(pack("i" * len(npart), *npart))
            f.write(pack("i" * len(fill), *fill))
            f.write(pack("i", nbytes))

            # Positions
            nbytes = int(len(pos) * 4)
            f.write(pack("i", nbytes))
            f.write(pack("f" * len(pos), *pos))
            f.write(pack("i", nbytes))

            # Velocities
            f.write(pack("i", nbytes))
            f.write(pack("f" * len(vel), *vel))
            f.write(pack("i", nbytes))

            # Ids
            nbytes = int(len(ids) * 4)
            f.write(pack("i", nbytes))
            f.write(pack("i" * len(ids), *ids))
            f.write(pack("i", nbytes))

            # Masses
            nbytes = len(mass) * 4
            f.write(pack("i", nbytes))
            f.write(pack("f" * len(mass), *mass))
            f.write(pack("i", nbytes))

            # Energy
            nbytes = len(u) * 4
            f.write(pack("i", nbytes))
            f.write(pack("f" * len(u), *u))
            f.write(pack("i", nbytes))

    elif format == 2:
        npart = array([ngas, 0, 0, 0, 0, 0])
        Nmass = array([0, 0, 0, 0, 0, 0])
        # Linearizing the 3D-array of the position and velocity
        pos = pos.ravel()
        vel = vel.ravel()

        dummy = zeros(npart[0])
        time = 0.0  # double
        redshift = 0.0  # double
        flag_sfr = 0  # int
        flag_feedback = 0  # int
        bytesleft = 256 - 6 * 4 - 6 * 8 - 8 - 8 - 2 * 4 - 6 * 4
        fill = zeros(int(bytesleft / 4.0), dtype=int)  # int

        with open(outfile, "wb") as f:

            nbytes = 256
            nbytes4 = 8

            f.write(pack("i", nbytes4))
            f.write(b"HEAD")
            f.write(pack("i", nbytes + 8))
            f.write(pack("i", nbytes4))
            # Header

            f.write(pack("i", nbytes))
            f.write(pack("i" * len(npart), *npart))
            f.write(pack("d" * len(Nmass), *Nmass))
            f.write(pack("d", time))
            f.write(pack("d", redshift))
            f.write(pack("i", flag_sfr))
            f.write(pack("i", flag_feedback))
            f.write(pack("i" * len(npart), *npart))
            f.write(pack("i" * len(fill), *fill))
            f.write(pack("i", nbytes))

            # Positions
            nbytes = int(len(pos) * 4)
            f.write(pack("i", nbytes4))
            f.write(b"POS ")
            f.write(pack("i", nbytes + 8))
            f.write(pack("i", nbytes4))

            f.write(pack("i", nbytes))
            f.write(pack("f" * len(pos), *pos))
            f.write(pack("i", nbytes))

            # Velocities
            f.write(pack("i", nbytes4))
            f.write(b"VEL ")
            f.write(pack("i", nbytes + 8))
            f.write(pack("i", nbytes4))

            f.write(pack("i", nbytes))
            f.write(pack("f" * len(vel), *vel))
            f.write(pack("i", nbytes))

            # Ids
            nbytes = int(len(ids) * 4)
            f.write(pack("i", nbytes4))
            f.write(b"ID  ")
            f.write(pack("i", nbytes + 8))
            f.write(pack("i", nbytes4))

            f.write(pack("i", nbytes))
            f.write(pack("i" * len(ids), *ids))
            f.write(pack("i", nbytes))

            # Masses
            nbytes = int(len(mass) * 4)
            f.write(pack("i", nbytes4))
            f.write(b"MASS")
            f.write(pack("i", nbytes + 8))
            f.write(pack("i", nbytes4))

            f.write(pack("i", nbytes))
            f.write(pack("f" * len(mass), *mass))
            f.write(pack("i", nbytes))

            # Energy
            nbytes = int(len(u) * 4)
            f.write(pack("i", nbytes4))
            f.write(b"U   ")
            f.write(pack("i", nbytes + 8))
            f.write(pack("i", nbytes4))

            f.write(pack("i", nbytes))
            f.write(pack("f" * len(u), *u))
            f.write(pack("i", nbytes))

    elif format == 3:
        with File(outfile, "w") as f:
            f.create_group("Header")
            f.create_group("PartType0")
            f["Header"].attrs["NumPart_ThisFile"] = array(
                [ngas, 0, 0, 0, 0, 0], dtype=uint32
            )
            f["Header"].attrs["NumPart_Total"] = array(
                [ngas, 0, 0, 0, 0, 0], dtype=uint32
            )
            f["Header"].attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
            f["Header"].attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            f["Header"].attrs["NumFilesPerSnapshot"] = 1
            f["Header"].attrs["Time"] = 0.0
            f["Header"].attrs["Redshift"] = 0.0
            f["Header"].attrs["BoxSize"] = boxsize
            f["Header"].attrs["Flag_Sfr"] = int32(0)
            f["Header"].attrs["Flag_Feedback"] = int32(0)
            f["Header"].attrs["Flag_Entropy_ICs"] = int32(0)

            f["PartType0"].create_dataset("Masses", data=mass.astype(float32))
            f["PartType0"].create_dataset("Coordinates", data=pos.astype(float32))
            f["PartType0"].create_dataset("Velocities", data=vel.astype(float32))
            f["PartType0"].create_dataset("ParticleIDs", data=ids.astype(int32))
            f["PartType0"].create_dataset("InternalEnergy", data=u.astype(float32))
            f["PartType0"].create_dataset("SmoothingLength", data=h.astype(float32))

    else:
        print("Format {} unknown or not implemented. Exiting.".format(format))
        exit()
