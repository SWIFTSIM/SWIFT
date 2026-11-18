#!/bin/bash
#
# Cosmological units and density definition of the cosmic-ray source fire2.
#
# Needs a SWIFT built with --with-cooling=geajay_4 and the program
# testGeajayKeys built from tests/testGeajayKeys.c. Not part of "make check". Usage:
#
#   tests/testGeajayCosmicRaysFire2Cosmology.sh <testGeajayKeys> <fixture dir> <work dir>
#
# The fixture dir holds CloudyData_UVB=HM2012.h5. The program passes one gas
# particle through the module with GEAJAY_KEYS_SCALE_FACTOR set, so that the
# particle density is comoving as in a cosmological run, and prints the
# cosmic-ray value that reaches the library (PARTCALL lines). The expected
# value is computed here from the command-line inputs only:
#   n_H = X (1 - Z) rho_comoving a^-3 / m_p,  value = e_cr n_H / (0.01 + n_H).
# Exit code 0: every check passes.

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

mk() {
  local IFS=';' l
  cp "$here/testGeajayKeys.yml" "$1"
  for l in $2; do [ -n "$l" ] && echo "  $l" >> "$1"; done
}

verdict() {
  checks=$((checks + 1))
  if [ "$2" -eq 1 ]; then
    printf '%s\tPASS\t%s\n' "$1" "$3"
  else
    printf '%s\tFAIL\t%s\n' "$1" "$3"
    fails=$((fails + 1))
  fi
}

ECR=2.5
CR="cosmic_rays: reduced;cosmic_ray_rate_per_energy_density: 1e-16;cosmic_ray_coulomb_loss_coefficient: 1e-16;cosmic_ray_hadronic_loss_coefficient: 1e-16;cosmic_ray_input: energy_density;cosmic_ray_source: fire2;cosmic_ray_fire2_energy_density: $ECR"

# run <name> <scale factor> <comoving density> : the particle at that a
run() {
  mk "$1.yml" "$CR"
  GEAJAY_KEYS_SCALE_FACTOR=$2 "$bin" "$1.yml" --particle "$3" 1e3 > "$1.out" 2>&1
  echo $? > "$1.status"
}

# value <name> <call> : the cosmic-ray value of a PARTCALL
value() { sed -n "s/^PARTCALL $2 cosmic_ray_value=\([^ ]*\) .*/\1/p" "$1.out" | head -1; }

# expect <a> <rho_comoving> <name> <definition>: the value of a definition
# definition: phys (X(1-Z) rho a^-3), comoving (X(1-Z) rho), mass (rho a^-3 / m_p)
expect() {
  local xh mp du z
  xh=$(sed -n 's/^MOD HydrogenFractionByMass=//p' "$3.out")
  mp=$(sed -n 's/^MOD proton_mass_cgs=//p' "$3.out")
  du=$(sed -n 's/^MOD density_units=//p' "$3.out")
  z=$(sed -n 's/^PARTCALL solve .* density=\([^ ]*\) metal=\(.*\)$/\1 \2/p' "$3.out" | head -1)
  awk -v a="$1" -v r="$2" -v e=$ECR -v x="$xh" -v mp="$mp" -v du="$du" -v def="$4" -v zz="$z" '
    BEGIN { split(zz, q, " "); zf = q[2] / q[1]
            if (def == "phys") n = x * (1 - zf) * r / (a * a * a) * du / mp
            else if (def == "comoving") n = x * (1 - zf) * r * du / mp
            else n = r / (a * a * a) * du / mp
            printf "%.17g", e * n / (0.01 + n) }'
}

# close <x> <y> <tolerance>: 1 when equal to the relative tolerance
close_to() { awk -v x="$1" -v y="$2" -v t="$3" 'BEGIN { r = (x - y) / y; exit !(r < t && r > -t) }' && echo 1 || echo 0; }

check_case() {
  local name=$1 a=$2 rho=$3 ok=1 c v w
  run "$name" "$a" "$rho"
  for c in cooling_time energy_rates solve; do
    v=$(value "$name" $c)
    w=$(expect "$a" "$rho" "$name" phys)
    [ -n "$v" ] && [ "$(close_to "$v" "$w" 1e-12)" -eq 1 ] || ok=0
  done
  verdict "fire2_cosmo_${name}_value_from_physical_n_H" $ok \
    "a=$a comoving rho=$rho: value $(value "$name" solve), expected $(expect "$a" "$rho" "$name" phys)"
}

# n_H near the 0.01 cm^-3 knee, where the factor is sensitive to n_H
check_case a05_knee 0.5 4e4
check_case a025_knee 0.25 5e3
# n_H far above the knee: the value tends to e_cr
check_case a05_dense 0.5 4e8

# the same physical gas gives the same value at any a
run a1_same_physical 1.0 3.2e5
verdict fire2_cosmo_value_depends_on_physical_density_only \
  "$(close_to "$(value a1_same_physical solve)" "$(value a05_knee solve)" 1e-12)" \
  "a=1 rho=3.2e5: $(value a1_same_physical solve); a=0.5 rho=4e4: $(value a05_knee solve)"

# the other density definitions give other values: the test would fail on them
obs=$(value a05_knee solve)
for def in comoving mass; do
  alt=$(expect 0.5 4e4 a05_knee $def)
  verdict "fire2_cosmo_definition_${def}_is_distinguished" \
    "$(awk -v o="$obs" -v w="$alt" 'BEGIN { r = (o - w) / w; print (r > 1e-2 || r < -1e-2) ? 1 : 0 }')" \
    "observed $obs, definition $def would give $alt"
done

# the value is passed in eV cm^-3, in no other unit: far above the knee it is e_cr
dense=$(value a05_dense solve)
verdict fire2_cosmo_value_tends_to_e_cr_in_ev_per_cm3 \
  "$(close_to "$dense" "$ECR" 1e-3)" "dense value $dense for e_cr $ECR eV cm^-3"

echo "checks: $checks failures: $fails"
[ $fails -eq 0 ]
