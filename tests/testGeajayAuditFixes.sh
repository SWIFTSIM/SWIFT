#!/bin/bash
#
# Start-up and report fixes of the geajay cooling module after the weakness
# audit: the count of the dust warnings of the library, the helium of the
# default start and the text of the restart refusals.
#
# Needs a SWIFT built with --with-cooling=geajay_3 and the program
# testGeajayKeys built from tests/testGeajayKeys.c. Not part of "make check".
# Usage:
#
#   tests/testGeajayAuditFixes.sh <testGeajayKeys> <fixture dir> <work dir>
#
# The fixture dir holds CloudyData_UVB=HM2012.h5. Exit code 0: every check
# passes.

here=$(dirname "$(realpath "$0")")
bin=$(realpath "$1")
fixture=$(realpath "$2")
mkdir -p "$3" || exit 2
out=$(realpath "$3")
unset LD_LIBRARY_PATH
cd "$out" || exit 2
ln -sf "$fixture/CloudyData_UVB=HM2012.h5" .

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

# The solves of a step with a library warning in their status are counted and
# reported by the next cooling_update. GEAJAY_KEYS_FORCE_FLAGS adds the bits.
flag_report() {
  local name=$1 bits=$2 dust=$3 metal=$4 line
  if [ -n "$bits" ]; then
    GEAJAY_KEYS_FORCE_FLAGS=$bits "$bin" "$here/testGeajayKeys.yml" --particle 1e-3 100 > "$name.out" 2>&1
  else
    "$bin" "$here/testGeajayKeys.yml" --particle 1e-3 100 > "$name.out" 2>&1
  fi
  line=$(grep -m1 "with a warning in the status" "$name.out")
  if [ "$dust$metal" = "00" ]; then
    [ -z "$line" ] && grep -q KEYPASS_DONE "$name.out"
    verdict "$name" $((1 - $?)) "no report line, as expected with no warning"
  else
    grep -q "status: $dust with the dust temperature not converged (the best estimate was used), $metal with the metal density" <<< "$(tr '\n' ' ' <<< "$line")"
    verdict "$name" $((1 - $?)) "$(sed 's/.*status: //' <<< "$line")"
  fi
}

flag_report dust_flags_none "" 0 0
flag_report dust_flag_not_converged 2 1 0
flag_report dust_flag_metal_above_density 1 0 1
flag_report dust_flags_both 3 1 1

# The restart refusals do not advise a start from a snapshot, which resets the species.
restore_case() {
  local kind=$1 st
  GEAJAY_KEYS_RESTORE=$kind "$bin" "$here/testGeajayKeys.yml" > "restore_$kind.out" 2>&1
  st=$?
  if [ "$kind" = round_trip ]; then
    [ $st -eq 0 ] && grep -q "RESTORE round_trip ok" "restore_$kind.out"
    verdict "restart_$kind" $((1 - $?)) "a valid restart file is read back"
  else
    [ $st -ne 0 ] && grep -qF "real physics never restarts from a snapshot" "restore_$kind.out" &&
      ! grep -qF "Start again from a snapshot" "restore_$kind.out"
    verdict "restart_$kind" $((1 - $?)) "exit $st: $(grep -o 'Start again from [a-z ]*' "restore_$kind.out" | head -n 1)"
  fi
}

restore_case round_trip
restore_case no_record
restore_case other_layout

echo "failures: $fails"
[ $fails -eq 0 ]
