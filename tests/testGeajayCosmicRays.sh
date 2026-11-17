#!/bin/bash
#
# Source modes of the cosmic-ray value of the geajay cooling module.
#
# Needs a SWIFT built with --with-cooling=geajay_3 or geajay_4 and the program
# testGeajayKeys built from tests/testGeajayKeys.c. Not part of "make check". Usage:
#
#   tests/testGeajayCosmicRays.sh <testGeajayKeys> <fixture dir> <work dir>
#
# The fixture dir holds CloudyData_UVB=HM2012.h5. The program reads a parameter
# file, then passes one gas particle through the module (cooling time, energy
# rates and one solve) and prints the cosmic-ray value that reaches the library
# at each call (PARTCALL lines). Exit code 0: every check passes.

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

# mk <file> <lines separated by ;>: the key harness file plus the lines
mk() {
  local IFS=';' l
  cp "$here/testGeajayKeys.yml" "$1"
  for l in $2; do [ -n "$l" ] && echo "  $l" >> "$1"; done
}

# verdict <name> <ok: 0 or 1> <evidence>
verdict() {
  checks=$((checks + 1))
  if [ "$2" -eq 1 ]; then
    printf '%s\tPASS\t%s\n' "$1" "$3"
  else
    printf '%s\tFAIL\t%s\n' "$1" "$3"
    fails=$((fails + 1))
  fi
}

# values <file> <call>: the cosmic_ray_value of every PARTCALL of that call
values() { sed -n "s/^PARTCALL $2 cosmic_ray_value=\([^ ]*\) .*/\1/p" "$1"; }

# all_equal <file> <expected>: the three calls pass the expected number
all_equal() {
  local c n=0 v
  for c in cooling_time energy_rates solve; do
    v=$(values "$1" $c)
    [ -n "$v" ] && awk -v a="$v" -v b="$2" 'BEGIN { exit !(a == b) }' && n=$((n + 1))
  done
  [ $n -eq 3 ]
}

CRB="cosmic_rays: reduced;cosmic_ray_rate_per_energy_density: 1e-16;cosmic_ray_coulomb_loss_coefficient: 1e-16;cosmic_ray_hadronic_loss_coefficient: 1e-16"
CR="$CRB;cosmic_ray_input: energy_density"

# run <name> <lines> <rho> <u> [value to set]
run() {
  mk "$1.yml" "$2"
  "$bin" "$1.yml" --particle "${@:3}" > "$1.out" 2>&1
  echo $? > "$1.status"
}
# refuse <name> <lines> <text>: the file must stop with the text
refuse() {
  mk "$1.yml" "$2"
  "$bin" "$1.yml" > "$1.out" 2>&1
  local st=$?
  verdict "$1" $([ $st -ne 0 ] && grep -qF -- "$3" "$1.out" && echo 1 || echo 0) \
    "refused: $3"
}

# constant: the value of the file reaches the library, nothing is stored
run constant "$CR;cosmic_ray_initial_value: 0.25" 1e-3 1e3
verdict constant_value_reaches_library \
  $(all_equal constant.out 0.25 && echo 1 || echo 0) "0.25 at the three calls"
verdict constant_stores_nothing \
  $(grep -qx "PART stored_value_after_first_init=0" constant.out && echo 1 || echo 0) \
  "stored value 0"
verdict constant_source_is_default \
  $(grep -qx "MOD cosmic_ray_source=0" constant.out && echo 1 || echo 0) "source 0"

# snapshot fields: the value of the run in constant, the stored value in per_particle
run snap_const "$CR;cosmic_ray_initial_value: 0.25" 1e-3 1e3
run snap_pp "$CR;cosmic_ray_source: per_particle" 1e-3 1e3 0.25
fc=$(sed -n 's/^FIELD after_first_init CosmicRayValues \(.*\)$/\1/p' snap_const.out)
fp=$(sed -n 's/^FIELD after_set CosmicRayValues \(.*\)$/\1/p' snap_pp.out)
verdict snapshot_value_is_the_run_value_in_constant \
  $([ -n "$fc" ] && [ "$fc" = "$fp" ] && ! grep -q "value=0$" <<<"$fc" && echo 1 || echo 0) \
  "CosmicRayValues: constant 0.25 and stored 0.25 give '$fc'"
verdict snapshot_fields_per_source \
  $(grep -q "^FIELD after_first_init CosmicRayEnergyLossRates" snap_pp.out && echo 1 || echo 0) \
  "the loss rate field is written"

# per_particle: starts at 0, then the value written through the function
run pp_start "$CR;cosmic_ray_source: per_particle" 1e-3 1e3
verdict per_particle_starts_at_zero \
  $(all_equal pp_start.out 0 && grep -qx "PART stored_value_after_first_init=0" pp_start.out && echo 1 || echo 0) \
  "0 at the three calls"
run pp_set "$CR;cosmic_ray_source: per_particle" 1e-3 1e3 3.5
verdict per_particle_written_value_reaches_library \
  $(all_equal pp_set.out 3.5 && echo 1 || echo 0) "3.5 at the three calls"
verdict per_particle_not_depleted_by_cooling \
  $(grep -qx "PART stored_value_after_solve=3.5" pp_set.out && echo 1 || echo 0) \
  "3.5 after the solve"
lg=$(sed -n 's/^PART loss_rate_getter=//p' pp_set.out)
lf=$(sed -n 's/^PART loss_rate_field=//p' pp_set.out)
verdict per_particle_loss_rate_readable \
  $(awk -v a="$lg" -v b="$lf" 'BEGIN { exit !(a > 0 && a == b) }' && echo 1 || echo 0) \
  "getter $lg, field $lf"
mk pp_ref.yml "$CR;cosmic_ray_initial_value: 1e-3"
"$bin" pp_ref.yml | grep "^LIB cr\." > pp_ref.lib
grep "^LIB cr\." pp_start.out > pp_start.lib
verdict per_particle_input_is_energy_density \
  $(cmp -s pp_ref.lib pp_start.lib && echo 1 || echo 0) \
  "LIB cr.* lines equal those of cosmic_ray_input: energy_density"
mk pp_noin.yml "cosmic_rays: reduced;cosmic_ray_source: per_particle;cosmic_ray_rate_per_energy_density: 1e-16;cosmic_ray_coulomb_loss_coefficient: 1e-16;cosmic_ray_hadronic_loss_coefficient: 1e-16"
"$bin" pp_noin.yml | grep "^LIB cr\." > pp_noin.lib
verdict per_particle_input_implied \
  $(cmp -s pp_ref.lib pp_noin.lib && echo 1 || echo 0) \
  "the input is energy_density without the key"

# fire2: e_cr n_H / (0.01 + n_H), from the density of the particle
ecr=2.5
for rho in 3.3e5 3.3e7; do
  run fire2_$rho "$CR;cosmic_ray_source: fire2;cosmic_ray_fire2_energy_density: $ecr" $rho 1e3
  xh=$(sed -n 's/^MOD HydrogenFractionByMass=//p' fire2_$rho.out)
  mp=$(sed -n 's/^MOD proton_mass_cgs=//p' fire2_$rho.out)
  du=$(sed -n 's/^MOD density_units=//p' fire2_$rho.out)
  ok=1
  for c in cooling_time energy_rates solve; do
    line=$(grep -m1 "^PARTCALL $c " fire2_$rho.out)
    d=$(sed -n 's/.* density=\([^ ]*\) .*/\1/p' <<<"$line")
    m=$(sed -n 's/.* metal=\(.*\)$/\1/p' <<<"$line")
    v=$(sed -n 's/.* cosmic_ray_value=\([^ ]*\) .*/\1/p' <<<"$line")
    awk -v d="$d" -v m="$m" -v v="$v" -v e=$ecr -v x="$xh" -v mp="$mp" -v du="$du" '
      BEGIN { nh = x * (d - m) * du / mp; w = e * nh / (0.01 + nh)
              r = (v - w) / w; exit !(r < 1e-12 && r > -1e-12) }' || ok=0
  done
  nh=$(awk -v d="$d" -v m="$m" -v x="$xh" -v mp="$mp" -v du="$du" 'BEGIN { printf "%.6g", x * (d - m) * du / mp }')
  verdict "fire2_value_at_n_H_$nh" $ok "e_cr $ecr eV cm^-3, value $v"
done
verdict fire2_snapshot_has_no_value_field \
  $(! grep -q "CosmicRayValues" fire2_3.3e7.out && grep -q "^FIELD after_first_init CosmicRayEnergyLossRates" fire2_3.3e7.out && echo 1 || echo 0) \
  "no CosmicRayValues, the loss rate is written"
verdict fire2_stores_nothing \
  $(grep -qx "PART stored_value_after_first_init=0" fire2_3.3e7.out && grep -qx "PART stored_value_after_solve=0" fire2_3.3e7.out && echo 1 || echo 0) \
  "stored value 0"
run fire2_set "$CR;cosmic_ray_source: fire2;cosmic_ray_fire2_energy_density: $ecr" 3.3e7 1e3 7
verdict fire2_refuses_a_written_value \
  $([ "$(cat fire2_set.status)" -ne 0 ] && grep -qF "cosmic_ray_source is not per_particle" fire2_set.out && echo 1 || echo 0) \
  "the function stops the run"

# keys
refuse source_unknown "$CR;cosmic_ray_source: columns" "must be constant, per_particle or fire2"
refuse fire2_needs_energy "$CR;cosmic_ray_source: fire2" "cosmic_ray_fire2_energy_density is required"
refuse fire2_energy_not_positive "$CR;cosmic_ray_source: fire2;cosmic_ray_fire2_energy_density: 0" "must be finite and > 0"
refuse energy_needs_fire2 "$CR;cosmic_ray_fire2_energy_density: 1" "cosmic_ray_source is not fire2"
refuse initial_value_with_per_particle "$CR;cosmic_ray_source: per_particle;cosmic_ray_initial_value: 1e-3" "cosmic_ray_source is not constant"
refuse initial_value_with_fire2 "$CR;cosmic_ray_source: fire2;cosmic_ray_fire2_energy_density: 1;cosmic_ray_initial_value: 1e-3" "cosmic_ray_source is not constant"
refuse per_particle_rate_input "$CRB;cosmic_ray_source: per_particle;cosmic_ray_input: ionisation_rate" "must be energy_density"
refuse source_without_cosmic_rays "cosmic_ray_source: per_particle" "cosmic_rays is off"

echo "checks: $checks failures: $fails"
[ $fails -eq 0 ]
