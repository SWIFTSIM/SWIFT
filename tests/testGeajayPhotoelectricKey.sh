#!/bin/bash
#
# Parse test of the GeajayCooling dust and photoelectric keys.
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
metallicity=""

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
  # metallicity, when set, replaces GEARChemistry:initial_metallicity
  [ -n "$metallicity" ] &&
    sed -i "s/initial_metallicity: .*/initial_metallicity: $metallicity/" "$dir"/coolingBox.yml
  (cd "$dir" && "$swift" --hydro --cooling --threads=1 -n 2 coolingBox.yml > out.log 2>&1)
  status=$?
}

# expect_ok <name> <chemistry_data values, as "member=value ..."> <extra lines>
expect_ok() {
  local name=$1 wanted=$2 item bad=""
  run "$name" "${@:3}"
  [ $status -ne 0 ] && bad="exit $status"
  for item in $wanted; do
    grep -q "chemistry_data.${item%%=*} = ${item#*=}\$" "$work/$name/out.log" ||
      bad="$bad, wanted ${item%%=*} = ${item#*=}"
  done
  if [ -n "$bad" ]; then
    echo "FAIL $name: $bad"
    failures=$((failures + 1))
  else
    echo "ok   $name ($wanted)"
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

expect_ok default photoelectric_heating=0 ""
expect_ok default_dust photoelectric_heating=2 "$basic"
expect_ok none_dust photoelectric_heating=0 "$basic;photoelectric_heating_model: none"
expect_ok fixed_rate photoelectric_heating=1 "$basic;photoelectric_heating_model: fixed_rate"
expect_ok constant_efficiency photoelectric_heating=2 "$basic;photoelectric_heating_model: constant_efficiency"
expect_ok grain_charge_epsilon photoelectric_heating=3 "$basic;photoelectric_heating_model: grain_charge_epsilon"
expect_ok density_epsilon photoelectric_heating=4 "$basic;photoelectric_heating_model: density_epsilon"
expect_ok number_4 photoelectric_heating=4 "$basic;photoelectric_heating_model: 4"
expect_ok default_name photoelectric_heating=2 "$basic;photoelectric_heating_model: default"
expect_ok no_dust_none photoelectric_heating=0 "photoelectric_heating_model: none"

expect_refused bad_name "photoelectric_heating_model" "$basic;photoelectric_heating_model: bogus"
expect_refused bad_high "must be default (-1), none (0)" "$basic;photoelectric_heating_model: 5"
expect_refused bad_low "must be default (-1), none (0)" "$basic;photoelectric_heating_model: -2"
expect_refused bad_huge "photoelectric_heating_model" "$basic;photoelectric_heating_model: 4294967297"

# SWIFT passes the metal mass fraction of every particle, so a dust switch
# alone is enough: the checks are the ones of the library.
expect_ok pe_no_dust photoelectric_heating=2 "photoelectric_heating_model: constant_efficiency"
expect_ok pe_fixed_no_dust photoelectric_heating=1 "photoelectric_heating_model: fixed_rate"
expect_ok pe_density_no_dust photoelectric_heating=4 "photoelectric_heating_model: density_epsilon"
expect_ok h2dust_network_1 h2_on_dust=1 "primordial_chemistry_network: H_He;h2_formation_on_dust: enabled"
expect_ok h2dust_network_0 h2_on_dust=1 "primordial_chemistry_network: tabulated;h2_formation_on_dust: enabled"
expect_refused h2dust_network_1_basic "h2_formation_on_dust is on with the network H_He" "primordial_chemistry_network: H_He;$basic;h2_formation_on_dust: enabled"
expect_ok remy_with_pe dust_to_gas_model=1 "dust_to_gas_model: remy_ruyer_2014;photoelectric_heating_model: constant_efficiency"
expect_ok remy_with_h2dust dust_to_gas_model=1 "dust_to_gas_model: remy_ruyer_2014;h2_formation_on_dust: enabled"
expect_ok remy_with_dust dust_to_gas_model=1 "dust_to_gas_model: remy_ruyer_2014;$basic"
expect_refused remy_no_switch "dust_to_gas_model = 1" "dust_to_gas_model: remy_ruyer_2014"

expect_ok rate_default photoelectric_heating_rate=8.5e-26 ""
expect_ok rate_minus_1 photoelectric_heating_rate=8.5e-26 "photoelectric_heating_model: fixed_rate;photoelectric_heating_rate: -1"
expect_ok rate_double photoelectric_heating_rate=1.7e-25 "photoelectric_heating_model: fixed_rate;photoelectric_heating_rate: 1.7e-25"
expect_ok rate_zero "photoelectric_heating=1 photoelectric_heating_rate=0" "photoelectric_heating_model: fixed_rate;photoelectric_heating_rate: 0"
expect_ok rate_unread_model photoelectric_heating_rate=3e-25 "photoelectric_heating_model: density_epsilon;photoelectric_heating_rate: 3e-25"
expect_refused rate_negative "photoelectric_heating_rate must be -1" "photoelectric_heating_model: fixed_rate;photoelectric_heating_rate: -2e-26"
expect_refused rate_inf "photoelectric_heating_rate must be -1" "photoelectric_heating_model: fixed_rate;photoelectric_heating_rate: inf"
expect_ok recomb_default dust_recombination_cooling=0 ""
expect_ok recomb_default_dust dust_recombination_cooling=1 "$basic"
expect_ok recomb_default_name dust_recombination_cooling=1 "$basic;dust_recombination_cooling: default"
expect_ok recomb_enabled dust_recombination_cooling=1 "dust_recombination_cooling: enabled"
expect_ok recomb_number dust_recombination_cooling=1 "dust_recombination_cooling: 1"
expect_ok recomb_disabled_dust dust_recombination_cooling=0 "$basic;dust_recombination_cooling: disabled"
expect_ok remy_with_recomb dust_to_gas_model=1 "dust_to_gas_model: remy_ruyer_2014;dust_recombination_cooling: enabled"
expect_refused recomb_bad_name "dust_recombination_cooling" "dust_recombination_cooling: bogus"
expect_refused recomb_bad_high "dust_recombination_cooling must be default (-1)" "dust_recombination_cooling: 2"
expect_refused recomb_bad_low "dust_recombination_cooling must be default (-1)" "dust_recombination_cooling: -2"
expect_ok ref_ratio_default local_dust_to_gas_ratio=0.00939 ""
expect_ok ref_ratio_minus_1 local_dust_to_gas_ratio=0.00939 "reference_dust_to_gas_ratio: -1"
expect_ok ref_ratio_set local_dust_to_gas_ratio=0.01 "reference_dust_to_gas_ratio: 0.01"
expect_refused ref_ratio_zero "reference_dust_to_gas_ratio must be -1" "reference_dust_to_gas_ratio: 0"
expect_refused ref_ratio_negative "reference_dust_to_gas_ratio must be -1" "reference_dust_to_gas_ratio: -0.5"
expect_refused ref_ratio_nan "reference_dust_to_gas_ratio must be -1" "reference_dust_to_gas_ratio: nan"
expect_refused ref_ratio_inf "reference_dust_to_gas_ratio must be -1" "reference_dust_to_gas_ratio: inf"
expect_refused rate_nan "photoelectric_heating_rate must be -1" "photoelectric_heating_model: fixed_rate;photoelectric_heating_rate: nan"

expect_refused from_field_name "SWIFT has no per-particle dust field" "dust_to_gas_model: from_field"
expect_refused from_field_values "linear (0) or remy_ruyer_2014 (1), got 2" "dust_to_gas_model: from_field"

# The dust column of the per_band attenuation needs a given H2 shielding length.
per_band="background_mode: per_band;isrf_background_model: constant;isrf_background_habing: 1;dust_background_attenuation: per_band;h2_shielding_length_method: given_length;isrf_from_field: 1;$basic"
expect_ok xs_default "local_dust_to_gas_ratio=0.00939" "$per_band"
expect_refused xs_lw_nan "dust_cross_section_lw_cm2 must be a finite number >= 0, got nan" "$per_band;dust_cross_section_lw_cm2: nan"
expect_refused xs_pe_nan "dust_cross_section_pe_cm2 must be a finite number >= 0, got nan" "$per_band;dust_cross_section_pe_cm2: nan"
expect_refused xs_lw_inf "dust_cross_section_lw_cm2 must be a finite number >= 0, got inf" "$per_band;dust_cross_section_lw_cm2: inf"
expect_refused xs_pe_negative "dust_cross_section_pe_cm2 must be a finite number >= 0, got -1e-21" "$per_band;dust_cross_section_pe_cm2: -1e-21"

# expect_heating <name> <positive|zero|nonzero>: the "Dust" component of
# CoolingRates in the first snapshot of the run <name>, for every particle.
# With no dust_model it is the rate of the one switch that is on, and positive
# means heating.
expect_heating() {
  local out
  out=$(python3 -I -c "
import sys, h5py, numpy as np
f = h5py.File(sys.argv[1], 'r')
names = [n.decode() for n in f['SubgridScheme/NamedColumns/CoolingRates'][:]]
rate = f['PartType0/CoolingRates'][:, names.index('Dust')]
test = {'positive': np.all(rate > 0), 'zero': np.all(rate == 0), 'nonzero': np.all(rate != 0)}
ok = bool(np.all(np.isfinite(rate)) and test[sys.argv[2]])
print(ok, rate.min(), rate.max())
" "$work/$1/coolingBox_0000.hdf5" "$2" 2>&1 | tail -n 1)
  if [ "${out%% *}" != "True" ]; then
    echo "FAIL $1: wanted Dust rate $2 in every particle, got: $out"
    failures=$((failures + 1))
  else
    echo "ok   $1 (Dust rate $2: ${out#True })"
  fi
}

# expect_ratio <name a> <name b> <ratio> [differ]: the "Dust" component of run b
# over that of run a, in every particle, to 1e-6 (the environment variable
# ratio_tol replaces it). With a 4th argument the ratio must differ from <ratio>
# by more than 1e-3 in every particle.
expect_ratio() {
  local out
  out=$(python3 -I -c "
import os, sys, h5py, numpy as np
def dust(path):
    f = h5py.File(path, 'r')
    names = [n.decode() for n in f['SubgridScheme/NamedColumns/CoolingRates'][:]]
    return f['PartType0/CoolingRates'][:, names.index('Dust')]
a, b = dust(sys.argv[1]), dust(sys.argv[2])
r = b / a
dev = np.abs(r / float(sys.argv[3]) - 1)
ok = np.all(dev < float(os.environ.get('ratio_tol', '1e-6'))) if len(sys.argv) < 5 else np.all(dev > 1e-3)
print(bool(ok), r.min(), r.max())
" "$work/$1/coolingBox_0000.hdf5" "$work/$2/coolingBox_0000.hdf5" "$3" $4 2>&1 | tail -n 1)
  if [ "${out%% *}" != "True" ]; then
    echo "FAIL $2 over $1: wanted Dust rate ratio $3, got: $out"
    failures=$((failures + 1))
  else
    echo "ok   $2 over $1 (Dust rate ratio $3: ${out#True })"
  fi
}

expect_heating no_dust_none zero
expect_heating pe_no_dust positive
expect_heating pe_fixed_no_dust positive
expect_heating pe_density_no_dust positive
expect_heating rate_zero zero
expect_ratio pe_fixed_no_dust rate_double 2
expect_heating recomb_default zero
expect_heating recomb_enabled nonzero

# The heating and the per_band attenuation use one dust-to-gas ratio. At
# Z = 0.002 (Z / Z_sun = 0.154, below the transition of remy_ruyer_2014) the
# model gives D = r f, with r the reference ratio and f = 0.0530. A linear run
# at Z = f Z_sun has the same D, so its Dust rate must equal that of the
# remy_ruyer_2014 run, with and without the PE attenuation. The large
# cross-section makes the attenuation matter: the thin and the thick rates
# differ.
zsun=0.01295
ref=0.00939
z_remy=0.002
z_linear=$(python3 -I -c "
import math, sys
zsun, z = float(sys.argv[1]), float(sys.argv[2])
x = z / zsun
print(repr(zsun * 10 ** (1.239 + 3.1 * math.log10(x))))
" $zsun $z_remy)
thick="dust_cross_section_pe_cm2: 3e-16"
thin="dust_cross_section_pe_cm2: 0"
remy="dust_to_gas_model: remy_ruyer_2014;reference_dust_to_gas_ratio: $ref"
linear="dust_to_gas_model: linear;reference_dust_to_gas_ratio: $ref"
metallicity=$z_remy
expect_ok att_remy_thin "" "$per_band;$thin;$remy"
expect_ok att_remy_thick "" "$per_band;$thick;$remy"
metallicity=$z_linear
expect_ok att_linear_thin "" "$per_band;$thin;$linear"
expect_ok att_linear_thick "" "$per_band;$thick;$linear"
metallicity=
# The runs differ in Z. Their thin rates differ by up to 1.3e-6 (measured, cause
# not isolated), against 0.84 for the thick pair of the old code.
export ratio_tol=1e-5
expect_ratio att_remy_thin att_linear_thin 1
expect_ratio att_remy_thick att_linear_thick 1
unset ratio_tol
expect_ratio att_remy_thin att_remy_thick 1 differ

echo "failures: $failures"
[ $failures -eq 0 ]
