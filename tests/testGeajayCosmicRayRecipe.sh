#!/bin/bash
#
# Depletion recipe of geajay_cooling.rst under -ffast-math.
#
# The recipe reads the loss rate of a particle with
# cooling_get_cosmic_ray_energy_loss_rate() and depletes the stored energy
# density. The field holds NaN after a failed solve and the getter returns 0.
# A finiteness test of the getter can be removed by -ffast-math, so only a build
# with the SWIFT flags shows that the getter works. Needs a SWIFT built with
# --with-cooling=geajay_3 or geajay_4 and the program testGeajayCosmicRayRecipe
# (make -C tests testGeajayCosmicRayRecipe GEAJAY_SRC=<geajay source>), which
# does not compile without -ffast-math. Not part of "make check". Usage:
#
#   tests/testGeajayCosmicRayRecipe.sh <testGeajayCosmicRayRecipe> <fixture dir> <work dir>
#
# The fixture dir holds CloudyData_UVB=HM2012.h5. Every check prints one line.
# Exit code 0: every check passes.

here=$(dirname "$(realpath "$0")")
bin=$(realpath "$1")
fixture=$(realpath "$2")
mkdir -p "$3" || exit 2
out=$(realpath "$3")
unset LD_LIBRARY_PATH
cd "$out" || exit 2
ln -sf "$fixture/CloudyData_UVB=HM2012.h5" .

cat > params.yml <<'YML'
InternalUnitSystem:
  UnitMass_in_cgs:     2.0e33
  UnitLength_in_cgs:   3.0857e21
  UnitVelocity_in_cgs: 1.0e5
  UnitCurrent_in_cgs:  1
  UnitTemp_in_cgs:     1
SPH:
  resolution_eta:        1.2348
  CFL_condition:         0.1
  minimal_temperature:   100.
GEARChemistry:
  initial_metallicity: 0.01295
GeajayCooling:
  dust_to_gas_model: linear
  cloudy_table_file: CloudyData_UVB=HM2012.h5
  with_uv_background: 0
  redshift: 0
  with_metal_cooling: 1
  thermal_time_myr: 5
  cosmic_rays: reduced
  cosmic_ray_input: energy_density
  cosmic_ray_rate_per_energy_density: 1e-16
  cosmic_ray_coulomb_loss_coefficient: 1e-16
  cosmic_ray_hadronic_loss_coefficient: 1e-16
  cosmic_ray_source: per_particle
YML

fails=0

# verdict <name> <ok: 0 or 1> <evidence>
verdict() {
  if [ "$2" -eq 1 ]; then
    printf '%s\tPASS\t%s\n' "$1" "$3"
  else
    printf '%s\tFAIL\t%s\n' "$1" "$3"
    fails=$((fails + 1))
  fi
}

# field <file> <name>: the value of before, getter, after or dt
field() { sed -n "s/^RECIPE .*\b$2=\([^ ]*\).*/\1/p" "$1"; }

# close <a> <b> <relative tolerance>
close() {
  awk -v a="$1" -v b="$2" -v t="$3" \
    'BEGIN { d = a - b; if (d < 0) d = -d; s = (b < 0 ? -b : b); exit !(d <= t * s) }'
}

# The step is long enough to deplete a real loss by about 12 percent.
DT=1e17

FAIL=1 "$bin" params.yml $DT > failed.out 2>&1
"$bin" params.yml $DT > real.out 2>&1
"$bin" params.yml 1e10 1e-11 > control.out 2>&1

cf=$(sed -n 's/^FASTMATH_TEST //p' failed.out)
case "$cf" in
  *-ffast-math*) verdict library_built_with_fast_math 1 "CFLAGS of libswiftsim hold -ffast-math" ;;
  *) verdict library_built_with_fast_math 0 "CFLAGS of libswiftsim: [$cf]" ;;
esac

b=$(field failed.out before); g=$(field failed.out getter); a=$(field failed.out after)
awk -v g="$g" 'BEGIN { exit !(g + 0 == 0 && g !~ /nan|inf/) }'
verdict failed_solve_getter_returns_zero $((1 - $?)) "getter=$g"
awk -v a="$a" -v b="$b" 'BEGIN { exit !(a != "" && a + 0 == b + 0 && b + 0 == 1) }'
verdict failed_solve_keeps_stored_energy $((1 - $?)) "stored $b before, $a after the recipe"

b=$(field real.out before); g=$(field real.out getter); a=$(field real.out after)
awk -v g="$g" 'BEGIN { exit !(g !~ /nan|inf/ && g + 0 > 0) }'
verdict real_solve_getter_is_finite_positive $((1 - $?)) "getter=$g"
exp=$(awk -v e="$b" -v l="$g" -v d="$DT" 'BEGIN { x = e - l * d; printf "%.9g", (x > 0 ? x : 0) }')
close "$a" "$exp" 1e-6 && ok=1 || ok=0
verdict real_solve_depletes_by_loss_times_step $ok "stored $a after the recipe, expected $exp"
awk -v a="$a" -v b="$b" 'BEGIN { exit !(a + 0 < b + 0 && a + 0 > 0) }'
verdict real_solve_depletion_is_visible $((1 - $?)) "stored $b before, $a after"

b=$(field control.out before); a=$(field control.out after)
close "$a" 0.9 1e-6 && ok=1 || ok=0
verdict control_loss_1e-11_depletes_to_0.9 $ok "stored $b before, $a after, expected 0.9"

echo "checks: 7 failures: $fails"
exit $((fails > 0))
