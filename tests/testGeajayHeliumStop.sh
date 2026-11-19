#!/bin/bash
#
# Test of the start-up refusal of helium ionisation rates combined with the HII
# coupling (GEARFeedback:HII_couple_ionization_rate) in the geajay cooling.
#
# Needs a SWIFT built with --with-cooling=geajay_X (X >= 1) and
# --with-subgrid-radiation=all or with the HII part, so it is not part of
# "make check". Usage:
#
#   tests/testGeajayHeliumStop.sh <swift binary> <fixture dir> <work dir>
#
# The fixture dir holds coolingBox.yml, coolingBox.hdf5, the Cloudy table of
# examples/Cooling/CoolingBox and the yields table radiation_fits_popII.hdf5.
# The yml needs a GeajayCooling block. Each case runs two steps with one thread. Exit code 0: every case behaves as listed.

swift=$(realpath "$1")
fixture=$(realpath "$2")
work=$3
mkdir -p "$work" || exit 2
work=$(realpath "$work")
unset LD_LIBRARY_PATH

failures=0
message="GEARFeedback:HII_couple_ionization_rate: 1 is combined with a nonzero"

# run <name> <HII_couple_ionization_rate> <extra GeajayCooling lines>: leaves
# the exit code in $status
run() {
  local name=$1 hii=$2
  shift 2
  local dir=$work/$name
  rm -rf "$dir"
  mkdir -p "$dir"
  cp "$fixture"/*.hdf5 "$fixture"/*.h5 "$dir"/
  awk -v extra="$*" '{print} /^GeajayCooling:/ {n = split(extra, a, ";"); for (i = 1; i <= n; i++) if (a[i] != "") print "  " a[i]}' \
    "$fixture"/coolingBox.yml > "$dir"/coolingBox.yml
  printf 'GEARFeedback:\n  HII_couple_ionization_rate: %s\n' "$hii" >> "$dir"/coolingBox.yml
  cat >> "$dir"/coolingBox.yml << 'EOF'
  yields_table: radiation_fits_popII.hdf5
  yields_table_first_stars: radiation_fits_popII.hdf5
  discrete_yields: 1
  elements: [Fe, Mg, O, C, Al, Ca, Ba, Zn, Eu]
  supernovae_Ia_energy_erg: 1e51
  supernovae_efficiency: 0
  with_stellar_wind_feedback: 0
EOF
  (cd "$dir" && "$swift" --hydro --cooling --threads=1 -n 2 coolingBox.yml > out.log 2>&1)
  status=$?
}

# expect_ok <name> <HII switch> <extra lines>: the run exits 0, the message is
# absent and the HII switch reached the module (the cooling banner shows it).
expect_ok() {
  run "$@"
  if [ $status -ne 0 ] || grep -q -- "$message" "$work/$1/out.log" ||
    ! grep -q "HII rate coupling $([ "$2" = 1 ] && echo on || echo off)," "$work/$1/out.log"; then
    echo "FAIL $1: exit $status, wanted a run with the HII coupling switch $2"
    failures=$((failures + 1))
  else
    echo "ok   $1 (runs)"
  fi
}

# expect_refused <name> <HII switch> <extra lines>: the run exits non-zero with
# the message, before the first step.
expect_refused() {
  run "$@"
  if [ $status -eq 0 ] || ! grep -q -- "$message" "$work/$1/out.log" ||
    grep -q "Starting the main loop\|\[00000\]" "$work/$1/out.log"; then
    echo "FAIL $1: exit $status, wanted a start-up refusal with \"$message\""
    failures=$((failures + 1))
  else
    echo "ok   $1 (refused, exit $status)"
  fi
}

fields="rt_use_fields: 1"
hei="rt_hei_ionisation_rate_per_s: 1e-17"
heii="rt_heii_ionisation_rate_per_s: 1e-17"

expect_ok hii_no_helium 1 "$fields;rt_hi_ionisation_rate_per_s: 1e-17"
expect_ok hii_zero_helium 1 "$fields;rt_hei_ionisation_rate_per_s: 0;rt_heii_ionisation_rate_per_s: 0"
expect_ok helium_no_hii 0 "$fields;$hei;$heii"
expect_refused hii_hei 1 "$fields;$hei"
expect_refused hii_heii 1 "$fields;$heii"
expect_refused hii_both 1 "$fields;$hei;$heii"

echo "failures: $failures"
[ $failures -eq 0 ]
