#!/bin/bash
#
# Run test of GeajayCooling:with_radiative_cooling on the uniform box of
# examples/Cooling/Network4Box (network H_He_H2_D, cold gas with cosmic rays).
# Needs a SWIFT built with --with-cooling=geajay_3 and python3 with h5py. Not
# part of "make check". Usage:
#
#   tests/testGeajayRadiativeCoolingKey.sh <swift binary> <Cloudy table> <work dir>
#
# Runs 100 steps with the key at 0, 1 and absent. With 0 the internal energy of
# every particle must equal its start value to 1e-12 (relative): the key
# removes every energy term of the library and nothing else changes the energy
# of a uniform static box. The bar is far below the float rounding (6e-8), so
# any energy update fails it. The H2 fraction must still grow by more than a
# factor 10. With 1 and absent the energy must grow by more than 10 percent
# (positive control) and the two runs must give the same energies and H2.
# Exit code 0: every check passes.

swift=$(realpath "$1")
table=$(realpath "$2")
mkdir -p "$3" || exit 2
work=$(realpath "$3")
box=$(realpath "$(dirname "$0")/../examples/Cooling/Network4Box")
unset LD_LIBRARY_PATH

cd "$work" || exit 2
cp "$box/makeIC.py" "$box/network3.yml" . || exit 2
ln -sf "$table" CloudyData_UVB=HM2012.h5
[ -e primordialBox.hdf5 ] || python3 makeIC.py > makeIC.log 2>&1 || exit 2

mk() {
  sed -e "s/basename: *network3/basename: $1/" \
    -e 's/time_end: .*/time_end: 2.08e-4/' \
    -e 's/delta_time: .*/delta_time: 2.08e-4/' \
    -e 's/^  dt_max:.*/  dt_max: 2.080078125e-6/' network3.yml > "$1.yml"
  [ -n "$2" ] && sed -i "/^GeajayCooling:/a\\  with_radiative_cooling: $2" "$1.yml"
  printf '\nGEARChemistry:\n  initial_metallicity: 0.\n' >> "$1.yml"
}
mk key0 0
mk key1 1
mk keyabs ""
for k in key0 key1 keyabs; do
  "$swift" --hydro --cooling --threads=1 $k.yml > $k.log 2>&1 ||
    { echo "FAIL $k: swift exit $?"; exit 1; }
done

python3 - <<'P'
import sys
import h5py
import numpy as np

def read(k, n):
    g = h5py.File(f"{k}_{n}.hdf5")["PartType0"]
    o = np.argsort(g["ParticleIDs"][:])
    return g["InternalEnergies"][:][o], g["H2I"][:][o]

u0, h0 = read("key0", "0000")
fails = 0

def check(name, ok, text):
    global fails
    print(("ok   " if ok else "FAIL ") + name + " (" + text + ")")
    fails += not ok

u, h = read("key0", "0001")
d = np.abs(u / u0 - 1).max()
check("key 0 keeps the energy", d <= 1e-12, f"max |du/u| = {d:.3g}, bar 1e-12")
r = h.mean() / h0.mean()
check("key 0 evolves H2", r > 10, f"H2 grows by {r:.4g}, bar 10")

u, h = read("key1", "0001")
g = u.mean() / u0.mean() - 1
check("key 1 changes the energy", g > 0.1, f"du/u = {g:.4g}, bar 0.1")
ua, ha = read("keyabs", "0001")
check("absent key equals key 1", np.array_equal(u, ua) and np.array_equal(h, ha),
      "energies and H2 identical")
sys.exit(1 if fails else 0)
P
