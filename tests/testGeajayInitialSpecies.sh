#!/bin/bash
#
# Initial species of the geajay cooling module: the H nuclei of the gas are
# shared between HI, HII, H-, H2, H2+ and HD, and a ratio set above the
# available hydrogen stops the run at start-up. The He nuclei are shared
# between HeI, HeII, HeIII and HeH+ in the same way when the HeI ratio is unset.
#
# Needs a SWIFT built with --with-cooling=geajay_3 and the program
# testGeajayKeys built from tests/testGeajayKeys.c. Not part of "make check".
# Usage:
#
#   tests/testGeajayInitialSpecies.sh <testGeajayKeys> <fixture dir> <work dir> [reference testGeajayKeys]
#
# The fixture dir holds CloudyData_UVB=HM2012.h5. The program starts one gas
# particle through the module and prints its species (SPECIES lines). With a
# reference program, built before the change, the species of a start without
# molecular ratio must be the same text. Exit code 0: every check passes.

here=$(dirname "$(realpath "$0")")
bin=$(realpath "$1")
fixture=$(realpath "$2")
mkdir -p "$3" || exit 2
out=$(realpath "$3")
ref=""
[ -n "$4" ] && ref=$(realpath "$4")
unset LD_LIBRARY_PATH
cd "$out" || exit 2
ln -sf "$fixture/CloudyData_UVB=HM2012.h5" .

fails=0

# mk <file> <lines separated by ;>: the key harness file plus the lines
mk() {
  local IFS=';' l
  command cp -f "$here/testGeajayKeys.yml" "$1"
  for l in $2; do [ -n "$l" ] && echo "  $l" >> "$1"; done
}

# verdict <name> <ok: 0 or 1> <evidence>
verdict() {
  if [ "$2" -eq 1 ]; then
    printf '%s\tPASS\t%s\n' "$1" "$3"
  else
    printf '%s\tFAIL\t%s\n' "$1" "$3"
    fails=$((fails + 1))
  fi
}

# run <program> <name> <lines>: write <name>.yml, run one particle, keep the SPECIES lines
run() {
  mk "$2.yml" "$3"
  "$1" "$2.yml" --particle 1e-3 100 > "$2.out" 2>&1
  grep '^SPECIES' "$2.out" | sed 's/^SPECIES [a-z_]* //' > "$2.species"
}

# sp <name> <species>: the value of one species
sp() { sed -n "s/^$2=//p" "$1.species"; }

# check_h <name> <HI> <H2I>: HI and H2I as expected, H nuclei sum is the hydrogen mass fraction
check_h() {
  local hi h2 sum
  hi=$(sp "$1" HI_frac)
  h2=$(sp "$1" H2I_frac)
  sum=$(awk -v hi="$hi" -v h2="$h2" -v hii="$(sp "$1" HII_frac)" -v hm="$(sp "$1" HM_frac)" \
    -v h2ii="$(sp "$1" H2II_frac)" -v hdi="$(sp "$1" HDI_frac)" \
    'BEGIN { printf "%.9g", hi + hii + hm + h2 + h2ii + hdi / 3 }')
  awk -v hi="$hi" -v h2="$h2" -v ehi="$2" -v eh2="$3" -v sum="$sum" -v nh="$NH" 'BEGIN {
    ok = (hi - ehi < 1e-6 * nh && ehi - hi < 1e-6 * nh && h2 - eh2 < 1e-6 * nh && eh2 - h2 < 1e-6 * nh &&
          sum - nh < 1e-6 * nh && nh - sum < 1e-6 * nh)
    exit !ok }'
  verdict "$1" $((1 - $?)) "HI=$hi (expect $2) H2I=$h2 (expect $3) H nuclei sum=$sum (expect $NH)"
}

# The hydrogen mass fraction that the particle gets without hydrogen_mass_fraction in the file
run "$bin" default ""
NH=$(sed -n 's/^HI_frac=//p' default.species)

run "$bin" h2_0p49 "initial_h2i_to_h_ratio: 0.49"
check_h h2_0p49 "$(awk -v nh="$NH" 'BEGIN { printf "%.9g", nh * (1 - 0.98) }')" \
  "$(awk -v nh="$NH" 'BEGIN { printf "%.9g", nh * 0.98 }')"

run "$bin" h2_1em4 "initial_h2i_to_h_ratio: 1e-4"
check_h h2_1em4 "$(awk -v nh="$NH" 'BEGIN { printf "%.9g", nh * (1 - 2e-4) }')" \
  "$(awk -v nh="$NH" 'BEGIN { printf "%.9g", nh * 2e-4 }')"

run "$bin" hii_h2 "initial_hii_to_h_ratio: 0.3;initial_h2i_to_h_ratio: 0.2"
check_h hii_h2 "$(awk -v nh="$NH" 'BEGIN { printf "%.9g", nh * (1 - 0.3 - 0.4) }')" \
  "$(awk -v nh="$NH" 'BEGIN { printf "%.9g", nh * 0.4 }')"

# Refusal: the ratios give more than the hydrogen of the gas
mk too_much.yml "initial_hii_to_h_ratio: 0.5;initial_h2i_to_h_ratio: 0.3"
"$bin" too_much.yml > too_much.out 2>&1
st=$?
grep -qF "hold more hydrogen than the gas has" too_much.out && [ $st -ne 0 ]
verdict too_much $((1 - $?)) "exit $st: $(grep -o 'which is [^ ]* > 1' too_much.out)"

# check_he <name> <HeI> <He total>: HeI as expected, He species sum to the helium of the gas
check_he() {
  local hei sum
  hei=$(sp "$1" HeI_frac)
  sum=$(awk -v a="$hei" -v b="$(sp "$1" HeII_frac)" -v c="$(sp "$1" HeIII_frac)" 'BEGIN { printf "%.9g", a + b + c }')
  awk -v hei="$hei" -v ehei="$2" -v sum="$sum" -v esum="$3" 'BEGIN {
    ok = (hei - ehei < 1e-6 * esum && ehei - hei < 1e-6 * esum && sum - esum < 1e-6 * esum && esum - sum < 1e-6 * esum)
    exit !ok }'
  verdict "$1" $((1 - $?)) "HeI=$hei (expect $2) He species sum=$sum (expect $3)"
}

# The helium mass fraction of the gas, with the hydrogen mass fraction of the default start
HE=$(awk -v nh="$NH" 'BEGIN { printf "%.9g", 1 - nh }')
HE_PER_H=$(awk -v nh="$NH" 'BEGIN { printf "%.9g", (1 - nh) / (4 * nh) }')

run "$bin" he_default ""
check_he he_default "$HE" "$HE"

run "$bin" he_heii "initial_heii_to_h_ratio: 0.01"
check_he he_heii "$(awk -v nh="$NH" -v he="$HE" 'BEGIN { printf "%.9g", he - 4 * nh * 0.01 }')" "$HE"

run "$bin" he_heii_heiii "initial_heii_to_h_ratio: 0.01;initial_heiii_to_h_ratio: 0.02"
check_he he_heii_heiii "$(awk -v nh="$NH" -v he="$HE" 'BEGIN { printf "%.9g", he - 4 * nh * 0.03 }')" "$HE"

# Refusal: HeII and HeIII hold more helium than the gas has
mk he_too_much.yml "initial_heii_to_h_ratio: 0.05;initial_heiii_to_h_ratio: 0.05"
"$bin" he_too_much.yml > he_too_much.out 2>&1
st=$?
grep -qF "hold more helium than the gas has" he_too_much.out && [ $st -ne 0 ]
verdict he_too_much $((1 - $?)) "exit $st: $(grep -o 'together hold [^ ]* He nuclei per H nucleus' he_too_much.out) (gas has $HE_PER_H)"

# HeH+ holds one He nucleus; the key exists with the network H_He_H2_D_ions only
run "$bin" he_hehii "initial_hehii_to_h_ratio: 0.01"
if [ -s he_hehii.species ]; then
  check_he he_hehii "$(awk -v nh="$NH" -v he="$HE" 'BEGIN { printf "%.9g", he - 4 * nh * 0.01 }')" \
    "$(awk -v nh="$NH" -v he="$HE" 'BEGIN { printf "%.9g", he - 4 * nh * 0.01 }')"
else
  printf 'he_hehii\tSKIP\tthis build has no HeH+ species\n'
fi

# An explicit HeI ratio is not reduced
run "$bin" he_explicit "initial_hei_to_h_ratio: 0.0689473687;initial_heii_to_h_ratio: 0.01"
check_he he_explicit "$(awk -v nh="$NH" 'BEGIN { printf "%.9g", 4 * nh * 0.0689473687 }')" "$HE"

# Without a helium ratio, or with an explicit HeI ratio, the species are those of the reference program
if [ -n "$ref" ]; then
  for c in "he_default|" "he_hei_explicit|initial_hei_to_h_ratio: 0.0789473687" \
    "he_hei_heii_explicit|initial_hei_to_h_ratio: 0.0689473687;initial_heii_to_h_ratio: 0.01"; do
    n=${c%%|*}
    run "$bin" "new_$n" "${c#*|}"
    run "$ref" "ref_$n" "${c#*|}"
    [ -s "ref_$n.species" ] && cmp -s "new_$n.species" "ref_$n.species"
    verdict "same_as_reference_$n" $((1 - $?)) "$(wc -l < "new_$n.species") species lines compared"
  done
  # A small HeII ratio: the reference program starts with the helium counted twice, below its 1e-3 sum check
  run "$ref" ref_he_small "initial_heii_to_h_ratio: 2e-4"
  run "$bin" new_he_small "initial_heii_to_h_ratio: 2e-4"
  [ -s ref_he_small.species ] && [ "$(sp ref_he_small HeI_frac)" = "$(sp he_default HeI_frac)" ] && [ "$(sp new_he_small HeI_frac)" != "$(sp he_default HeI_frac)" ]
  verdict reference_double_counts_helium $((1 - $?)) "HeI is $(sp ref_he_small HeI_frac) in the reference program and $(sp new_he_small HeI_frac) in the new one"
fi

# Without a molecular ratio the species are those of the reference program
if [ -n "$ref" ]; then
  for c in "default|" "hii_only|initial_hii_to_h_ratio: 0.001" "zero_h2|initial_h2i_to_h_ratio: 0"; do
    n=${c%%|*}
    run "$bin" "new_$n" "${c#*|}"
    run "$ref" "ref_$n" "${c#*|}"
    [ -s "ref_$n.species" ] && cmp -s "new_$n.species" "ref_$n.species"
    verdict "same_as_reference_$n" $((1 - $?)) "$(wc -l < "new_$n.species") species lines compared"
  done
fi

echo "failures: $fails"
[ $fails -eq 0 ]
