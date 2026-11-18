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
# cosmic-ray value that reaches the library (PARTCALL lines). The run uses
# "redshift: -1", so the library scale factor follows a as in a cosmological
# run. The expected value is computed here from the command-line inputs, the
# fixture metal fraction Z = (double)0.01f set by run_particle, and X, m_p and
# the density unit printed by the program (MOD lines):
#   n_H = X (1 - Z) rho_comoving a^-3 / m_p,  value = e_cr n_H / (0.01 + n_H).
# With GEAJAY_KEYS_MEAN_BARYON_DENSITY set (mean baryonic density today, g cm^-3),
# the FIRE-2 cut applies: the value is 0 unless n_H > 1000 rho_b(a) / m_p, with
# rho_b(a) = rho_b0 a^-3.
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
  sed 's/^\( *redshift:\).*/\1 -1/' "$here/testGeajayKeys.yml" > "$1"
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
# the metal mass fraction run_particle sets: (double)0.01f
ZFIX=0.009999999776482582
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
  local xh mp du
  xh=$(sed -n 's/^MOD HydrogenFractionByMass=//p' "$3.out")
  mp=$(sed -n 's/^MOD proton_mass_cgs=//p' "$3.out")
  du=$(sed -n 's/^MOD density_units=//p' "$3.out")
  awk -v a="$1" -v r="$2" -v e=$ECR -v x="$xh" -v mp="$mp" -v du="$du" -v def="$4" -v zf="$ZFIX" '
    BEGIN { if (a + 0 <= 0) { print "scale factor must be positive" > "/dev/stderr"; exit 2 }
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

# The FIRE-2 cut: value 0 unless n_H > 1000 x the mean baryonic density (cm^-3 of m_p)
RHOB0=4.2e-31

# run_cut <name> <scale factor> <comoving density>: as run, with a mean baryonic density
run_cut() {
  mk "$1.yml" "$CR"
  GEAJAY_KEYS_SCALE_FACTOR=$2 GEAJAY_KEYS_MEAN_BARYON_DENSITY=$RHOB0 "$bin" "$1.yml" --particle "$3" 1e3 > "$1.out" 2>&1
  echo $? > "$1.status"
}

# expect_cut <a> <rho_comoving> <name>: the value with the cut
expect_cut() {
  local xh mp du
  xh=$(sed -n 's/^MOD HydrogenFractionByMass=//p' "$3.out")
  mp=$(sed -n 's/^MOD proton_mass_cgs=//p' "$3.out")
  du=$(sed -n 's/^MOD density_units=//p' "$3.out")
  awk -v a="$1" -v r="$2" -v e=$ECR -v x="$xh" -v mp="$mp" -v du="$du" -v zf="$ZFIX" -v rb="$RHOB0" '
    BEGIN { n = x * (1 - zf) * r / (a * a * a) * du / mp
            cut = 1000 * rb / (a * a * a) / mp
            printf "%.17g", (n > cut) ? e * n / (0.01 + n) : 0 }'
}

# cut_case <name> <a> <rho_comoving> <zero|nonzero>
cut_case() {
  local name=$1 a=$2 rho=$3 want=$4 ok=1 c v w
  run_cut "$name" "$a" "$rho"
  w=$(expect_cut "$a" "$rho" "$name")
  for c in cooling_time energy_rates solve; do
    v=$(value "$name" $c)
    if [ "$want" = zero ]; then
      [ -n "$v" ] && [ "$(awk -v v="$v" 'BEGIN { print (v == 0) ? 1 : 0 }')" -eq 1 ] || ok=0
    else
      [ -n "$v" ] && [ "$(close_to "$v" "$w" 1e-12)" -eq 1 ] && [ "$(awk -v v="$v" 'BEGIN { print (v > 0) ? 1 : 0 }')" -eq 1 ] || ok=0
    fi
  done
  verdict "fire2_cut_${name}" $ok "a=$a comoving rho=$rho: value $(value "$name" solve), expected $w ($want)"
}

# n_H = 9.8e-3 at the knee, cut at a=0.5 is 2.0e-3 and at a=1 is 2.5e-4 (cm^-3)
cut_case a05_above_cut 0.5 4e4 nonzero
# 1e4 gives n_H = 2.4e-3, between the cut 2.0e-3 and the cut 2.6e-3 of a mean density divided by X
cut_case a05_just_above_cut 0.5 1e4 nonzero
cut_case a05_below_cut 0.5 6e3 zero
# the same physical gas (6e3 / a^3 = 4.8e4) is above the cut at a=1, where the mean density is 8 times lower
cut_case a1_same_physical_gas_above_cut 1.0 4.8e4 nonzero
cut_case a05_dense_above_cut 0.5 4e8 nonzero

echo "checks: $checks failures: $fails"
[ $fails -eq 0 ]
