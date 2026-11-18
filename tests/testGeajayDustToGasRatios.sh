#!/bin/bash
#
# Snapshot field DustToGasRatios of the geajay cooling module.
#
# Needs a SWIFT built with --with-cooling=geajay_3 or geajay_4 and the program
# testGeajayKeys built from tests/testGeajayKeys.c. Not part of "make check". Usage:
#
#   tests/testGeajayDustToGasRatios.sh <testGeajayKeys> <fixture dir> <work dir>
#
# The fixture dir holds CloudyData_UVB=HM2012.h5. The program passes one gas
# particle with a given metal mass fraction (GEAJAY_KEYS_METAL_FRACTION)
# through the module and prints the value of the snapshot field (FIELD lines).
# The expected value is computed here from the fit of Remy-Ruyer et al. (2014)
# and the linear law, from the printed local_dust_to_gas_ratio and solar metal
# fraction, without the library. Exit code 0: every check passes.

here=$(dirname "$(realpath "$0")")
bin=$(realpath "$1")
fixture=$(realpath "$2")
mkdir -p "$3" || exit 2
out=$(realpath "$3")
unset LD_LIBRARY_PATH
cd "$out" || exit 2
ln -sf "$fixture/CloudyData_UVB=HM2012.h5" .

fails=0
checks=0

# mk <file> <lines separated by ;>: the key harness file, a later line replaces
# the key of an earlier one
mk() {
  local IFS=';' l k
  cp "$here/testGeajayKeys.yml" "$1"
  for l in $2; do
    [ -z "$l" ] && continue
    k=${l%%:*}
    sed -i "/^  $k:/d" "$1"
    echo "  $l" >> "$1"
  done
}

verdict() {
  checks=$((checks + 1))
  printf '%s\t%s\t%s\n' "$1" "$([ "$2" -eq 1 ] && echo PASS || echo FAIL)" "$3"
  [ "$2" -eq 1 ] || fails=$((fails + 1))
}

# run <name> <lines> <metal mass fraction>
run() {
  mk "$1.yml" "$2"
  GEAJAY_KEYS_METAL_FRACTION=$3 "$bin" "$1.yml" --particle 1e-3 1e3 > "$1.out" 2>&1
}

# expected <law> <metal fraction> <fgr> <z_sun>: D of the law, double precision
expected() {
  awk -v law="$1" -v z="$2" -v fgr="$3" -v zs="$4" 'BEGIN {
    r = z / zs
    # Remy-Ruyer et al. (2014), broken power law of the gas-to-dust ratio:
    # log(G/D) = 2.21 - log(Z/Zsun) above the transition 12+log(O/H) = 8.10,
    # slope -1; below it the slope is -3.10. Solar value 8.69.
    a = 2.21; xt = 8.10 - 8.69; slope_low = 3.10
    d = fgr * r
    if (law == "rr" && r > 0) {
      lx = log(r) / log(10)
      if (lx <= xt)
        d = fgr * 10 ^ (xt - slope_low * (xt - lx))
    }
    printf "%.10g", d }'
}

field() { sed -n 's/^FIELD after_first_init DustToGasRatios .* value=//p' "$1"; }
zused() { sed -n 's/^PART metal_mass_fraction=//p' "$1"; }

# compare <name> <value> <expected>: relative difference within float storage
agree() {
  awk -v a="$1" -v b="$2" 'BEGIN { d = (a - b) / b; exit !(b > 0 && d < 1e-6 && d > -1e-6) }'
}

DUST="dust_model: basic"
run probe "$DUST;dust_to_gas_model: linear" 0.01
fgr=$(sed -n 's/^MOD local_dust_to_gas_ratio=//p' probe.out)
zs=$(sed -n 's/^MOD SolarMetalFractionByMass=//p' probe.out)

# 0.1 Zsun and 0.001 lie below the transition (0.257 Zsun), 0.01 above it
zlow=$(awk -v zs="$zs" 'BEGIN { printf "%.8g", 0.1 * zs }')
for z in $zlow 0.001 0.01; do
  for law in linear remy_ruyer_2014; do
    n=${law}_z$z
    run "$n" "$DUST;dust_to_gas_model: $law" "$z"
    zu=$(zused "$n.out")
    want=$(expected "$([ "$law" = linear ] && echo lin || echo rr)" "$zu" "$fgr" "$zs")
    got=$(field "$n.out")
    verdict "snapshot_dust_to_gas_${law}_at_metal_fraction_$z" \
      $([ -n "$got" ] && agree "$got" "$want" && echo 1 || echo 0) \
      "field $got, expected $want from Z $zu"
  done
done

# the two laws differ below the transition and agree above it
rl=$(field "linear_z$zlow.out")
rr=$(field "remy_ruyer_2014_z$zlow.out")
verdict ratio_of_laws_at_0.1_solar \
  $(awk -v a="$rr" -v b="$rl" 'BEGIN { r = a / b; exit !(r > 0.1 && r < 0.2) }' && echo 1 || echo 0) \
  "remy_ruyer_2014 / linear = $(awk -v a="$rr" -v b="$rl" 'BEGIN { printf "%.4g", a / b }')"
rl=$(field linear_z0.01.out)
rr=$(field remy_ruyer_2014_z0.01.out)
verdict laws_agree_above_transition \
  $([ "$rl" = "$rr" ] && echo 1 || echo 0) "linear $rl, remy_ruyer_2014 $rr"

# no metals: the ratio is 0
run zero_metals "$DUST;dust_to_gas_model: remy_ruyer_2014" 0
verdict zero_metals_give_zero_ratio \
  $([ "$(field zero_metals.out)" = 0 ] && echo 1 || echo 0) "field $(field zero_metals.out)"

# no dust switch: the field is not written (the library refuses remy_ruyer_2014 here)
run no_dust "dust_to_gas_model: linear" 0.001
verdict no_dust_switch_writes_no_field \
  $(! grep -q DustToGasRatios no_dust.out && grep -q KEYPASS_DONE no_dust.out && echo 1 || echo 0) \
  "no DustToGasRatios line"

# the per_band attenuation reads the ratio, so it also switches the field on
run per_band "background_mode: per_band;h2_shielding_length_method: given_length;isrf_background_model: constant;isrf_background_habing: 1;isrf_from_field: 1;dust_background_attenuation: per_band;dust_to_gas_model: linear" 0.001
verdict per_band_attenuation_writes_field \
  $([ -n "$(field per_band.out)" ] && echo 1 || echo 0) "field $(field per_band.out)"

echo "checks: $checks failures: $fails"
[ $fails -eq 0 ]
