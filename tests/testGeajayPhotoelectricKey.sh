#!/bin/bash
#
# Parse test of GeajayCooling:photoelectric_heating_model.
#
# Needs a SWIFT built with --with-cooling=geajay_X, so it is not part of
# "make check". Usage:
#
#   tests/testGeajayPhotoelectricKey.sh <swift binary> <fixture dir> <work dir>
#
# The fixture dir holds coolingBox.yml, coolingBox.hdf5 and the Cloudy table of
# examples/Cooling/CoolingBox. The yml needs a GeajayCooling block. Each case
# runs two steps with one thread. Exit code 0: every case behaves as listed.

swift=$(realpath "$1")
fixture=$(realpath "$2")
work=$3
mkdir -p "$work" || exit 2
work=$(realpath "$work")
unset LD_LIBRARY_PATH

failures=0

# run <name> <extra GeajayCooling lines>: leaves the exit code in $status
run() {
  local name=$1
  shift
  local dir=$work/$name
  rm -rf "$dir"
  mkdir -p "$dir"
  cp "$fixture"/coolingBox.hdf5 "$fixture"/*.h5 "$dir"/
  # The new lines go right under the GeajayCooling: header
  awk -v extra="$*" '{print} /^GeajayCooling:/ {n = split(extra, a, ";"); for (i = 1; i <= n; i++) print "  " a[i]}' \
    "$fixture"/coolingBox.yml > "$dir"/coolingBox.yml
  (cd "$dir" && "$swift" --hydro --cooling --threads=1 -n 2 coolingBox.yml > out.log 2>&1)
  status=$?
}

# expect_ok <name> <resolved value> <extra lines>
expect_ok() {
  run "$1" "${@:3}"
  if [ $status -ne 0 ] ||
    ! grep -q "chemistry_data.photoelectric_heating = $2\$" "$work/$1/out.log"; then
    echo "FAIL $1: exit $status, wanted photoelectric_heating = $2"
    failures=$((failures + 1))
  else
    echo "ok   $1 (photoelectric_heating = $2)"
  fi
}

# expect_refused <name> <text of the message> <extra lines>
expect_refused() {
  run "$1" "${@:3}"
  if [ $status -eq 0 ] || ! grep -q -- "$2" "$work/$1/out.log"; then
    echo "FAIL $1: exit $status, wanted a refusal with \"$2\""
    failures=$((failures + 1))
  else
    echo "ok   $1 (refused)"
  fi
}

basic="dust_model: basic"

expect_ok default 0 ""
expect_ok default_dust 2 "$basic"
expect_ok none_dust 0 "$basic;photoelectric_heating_model: none"
expect_ok fixed_rate 1 "$basic;photoelectric_heating_model: fixed_rate"
expect_ok constant_efficiency 2 "$basic;photoelectric_heating_model: constant_efficiency"
expect_ok grain_charge_epsilon 3 "$basic;photoelectric_heating_model: grain_charge_epsilon"
expect_ok density_epsilon 4 "$basic;photoelectric_heating_model: density_epsilon"
expect_ok number_4 4 "$basic;photoelectric_heating_model: 4"
expect_ok default_name 2 "$basic;photoelectric_heating_model: default"
expect_ok no_dust_none 0 "photoelectric_heating_model: none"

expect_refused bad_name "photoelectric_heating_model" "$basic;photoelectric_heating_model: bogus"
expect_refused bad_high "must be default (-1), none (0)" "$basic;photoelectric_heating_model: 5"
expect_refused bad_low "must be default (-1), none (0)" "$basic;photoelectric_heating_model: -2"
expect_refused bad_huge "photoelectric_heating_model" "$basic;photoelectric_heating_model: 4294967297"
expect_refused no_dust_model "dust_model is none" "photoelectric_heating_model: constant_efficiency"

echo "failures: $failures"
[ $failures -eq 0 ]
