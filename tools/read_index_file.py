#!/usr/bin/env python3

import sys
import numpy as np
import mmap
import struct
import os


class DarkMatterParticle():
    def __init__(self, pos, vel, acc, mass, ids):
        self.pos = pos
        self.vel = vel
        self.acc = acc
        self.mass = mass
        self.ids = ids

class LogfileReader():
    n_type = 7
    swift_type_dark_matter = 1

    # dtype for the particle's data
    index_dtype = np.dtype([("ids", np.ulonglong),
                            ("offset", np.uint64)])

    def __init__(self, basename):
        self.basename = basename
        logfile = basename + ".dump"
        self.f = open(logfile, "rb")
        self.m = mmap.mmap(self.f.fileno(), 0, prot=mmap.PROT_READ)
        self.read_logfile_header()

    def __del__(self):
        self.m.close()
        self.f.close()

    def read_int(self):
        val = self.m.read(4)
        return int.from_bytes(val, byteorder="little")

    def read_uint(self):
        val = self.m.read(4)
        return struct.unpack("I", val)[0]

    def read_offset(self):
        offset = self.m.read(6)
        return int.from_bytes(offset, byteorder="little")

    def read_string(self):
        string = self.m.read(self.string_length)
        string = struct.unpack("%is" % self.string_length, string)[0]
        i = string.find(b"\x00")
        return string[:i]

    def read_logfile_header(self):
        self.m.seek(0)
        self.vers = self.m.read(20)
        self.major = self.read_int()
        self.minor = self.read_int()
        self.offset_direction = self.read_int()
        self.first_record = self.read_offset()
        self.string_length = self.read_uint()
        self.n_masks = self.read_uint()

        mask_dtype = np.dtype([("name", "U%i" % self.string_length),
                               ("value", "int32"), ("size", "uint32")])

        self.masks = np.empty(self.n_masks, dtype=mask_dtype)
        for i in range(self.n_masks):
            self.masks[i]["name"] = self.read_string()
            self.masks[i]["value"] = 1 << i
            self.masks[i]["size"] = self.read_uint()

        print("Fields found (name, mask, size): ", self.masks)


    def read_record_header(self):
        mask = self.m.read(2)
        mask = struct.unpack("h", mask)[0]

        if (mask != 316):
            print(mask)

        offset = self.read_offset()
        return mask, offset

    def read_index(self, index_number):
        filename = basename + "_%04i.index" % index_number

        # Read the file
        with open(filename, "rb") as f:
            # read the time
            time = np.fromfile(f, dtype=float, count=1)
            time_int = np.fromfile(f, dtype=np.longlong, count=1)
            print("Time: {}, integer time: {}".format(
                time[0], time_int[0]))

            # read the number of particles
            nparts = np.fromfile(f, dtype=np.uint64, count=LogfileReader.n_type)

            print("Number of particles:", nparts)

            # read if the file is sorted
            sort = np.fromfile(f, dtype=np.bool, count=1)
            print("File is sorted?", sort[0])

            # read the memory alignment garbage
            n = ((f.tell() + 7) & ~7) - f.tell()
            f.read(n)

            # read the particles
            print("Particles data (ids / offset):")
            data = []
            for n in nparts:
                if n == 0:
                    data.append(None)
                    continue

                tmp = np.fromfile(f, dtype=LogfileReader.index_dtype, count=n)
                data.append(tmp)

            print("\t", data)

            # print the history of new particles
            n_new = np.fromfile(f, dtype=np.uint64, count=LogfileReader.n_type)
            print("New particles: ", n_new)

            new = []
            for n in n_new:
                if n == 0:
                    new.append(None)
                    continue

                new.append(np.fromfile(f, dtype=LogfileReader.index_dtype, count=n))
            print("\t", new)

            # print the history of particles removed
            n_rem = np.fromfile(f, dtype=np.uint64, count=LogfileReader.n_type)
            print("Particles removed: ", n_rem)

            removed = []
            for n in n_rem:
                if n == 0:
                    removed.append(None)
                    continue

                removed.append(np.fromfile(f, dtype=LogfileReader.index_dtype, count=n))
            print("\t", removed)
            return data, new, removed

    def seek(self, new_offset):
        self.m.seek(new_offset)

    def read_particle(self, mask):
        flag = None
        ref = self.masks["name"] == "SpecialFlags"
        if (mask & self.masks[ref]["value"]):
            flag = self.m.read(4)
            print("Flag: ", struct.unpack("i", flag)[0])

        pos = None
        ref = self.masks["name"] == "Coordinates"
        if (mask & self.masks[ref]["value"]):
            pos = self.m.read(24)
            pos = struct.unpack("ddd", pos)[0]

        vel = None
        ref = self.masks["name"] == "Velocities"
        if (mask & self.masks[ref]["value"]):
            vel = self.m.read(12)
            vel = struct.unpack("fff", vel)[0]

        acc = None
        ref = self.masks["name"] == "Accelerations"
        if (mask & self.masks[ref]["value"]):
            acc = self.m.read(12)
            acc = struct.unpack("fff", acc)[0]

        # print(pos, vel, acc)
        mass = None
        ref = self.masks["name"] == "Masses"
        if (mask & self.masks[ref]["value"]):
            mass = self.m.read(4)
            mass = struct.unpack("f", mass)[0]

        ids = None
        ref = self.masks["name"] == "ParticleIDs"
        if (mask & self.masks[ref]["value"]):
            ids = self.m.read(8)
            ids = struct.unpack("q", ids)[0]

        part = DarkMatterParticle(pos, vel, acc, mass, ids)
        return part


if __name__ == "__main__":
    logfile = sys.argv[-1]
    basename = logfile[:-5]

    # Open the logfile
    reader = LogfileReader(basename)

    # Read the first index file
    data, _, _ = reader.read_index(0)

    # Keep only the DM particles
    data = data[LogfileReader.swift_type_dark_matter]

    for i in range(len(data)):
        if i > 10:
            print("Early exit due to low performances")
            break

        off = data["offset"][i]
        off = int(off)
        index_ids = data["ids"][i]
        reader.seek(off)

        for i in range(10):
            mask, offset = reader.read_record_header()
            print(mask, offset, off)

            part = reader.read_particle(mask)
            if (index_ids != part.ids):
                raise Exception("Failed", index_ids, part.ids)

            if offset == 0:
                print("End of particle")
                break

            off += offset
            reader.seek(off)
