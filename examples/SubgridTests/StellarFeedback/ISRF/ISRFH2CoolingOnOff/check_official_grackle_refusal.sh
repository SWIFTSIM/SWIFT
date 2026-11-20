#!/bin/bash
# Check that a SWIFT built with the official Grackle (--with-cooling=grackle_2
# or grackle_3 and the ISRF compiled in) stops at start-up when
# GEARFeedback:ISRF_H2_line_absorption is 1, and starts with it 0.
#
# Usage: check_official_grackle_refusal.sh <swift> <ICs.hdf5> <cloudy.h5> <yields.hdf5> <workdir>
# Exit status: 0 if the switch-1 run exits non-zero with the refusal text and
# the switch-0 run does not print it and gets past start-up, 1 otherwise.

set -u
swift=$(realpath "$1"); ics=$(realpath "$2"); cloudy=$(realpath "$3")
yields=$(realpath "$4"); work=$5
here=$(dirname "$(realpath "$0")")
params="$here/../../../../SinkParticles/HomogeneousBox/params.yml"
refusal="ISRF_H2_line_absorption: 1 is not supported: Grackle shields the transported LW rate"
mkdir -p "$work"
cp -f "$ics" "$work/ICs_homogeneous_box.hdf5"
cp -f "$cloudy" "$yields" "$work/"

run_leg() {
  local sw=$1 d="$work/switch$1"
  mkdir -p "$d"
  cp -f "$ics" "$work/ICs_homogeneous_box.hdf5" "$d/" 2>/dev/null
  cp -f "$cloudy" "$yields" "$d/"
  sed -e "s#CloudyData_UVB=HM2012_high_density.h5#$(basename "$cloudy")#" \
      -e "s#^\( *with_interstellar_radiation_field:\) 0#\1 1#" \
      -e "s#^\( *H2_self_shielding:\) 0#\1 3#" \
      -e "s#^\( *ISRF_propagation:\) 0#\1 1#" \
      -e "/^ *with_interstellar_radiation_field:/a\\  ISRF_H2_line_absorption: $sw" \
      "$params" > "$d/params.yml"
  (cd "$d" && env -u LD_LIBRARY_PATH timeout 120 "$swift" --hydro --stars \
     --self-gravity --feedback --cooling --sync --limiter --threads=2 -n 1 \
     params.yml > output.log 2>&1)
  echo $? > "$d/rc"
}

run_leg 1
run_leg 0
rc1=$(cat "$work/switch1/rc"); rc0=$(cat "$work/switch0/rc")
n1=$(grep -c "$refusal" "$work/switch1/output.log")
n0=$(grep -c "$refusal" "$work/switch0/output.log")
echo "switch 1: exit $rc1, refusal text lines $n1"
echo "switch 0: exit $rc0, refusal text lines $n0, reached main: done $(grep -c "main: done" "$work/switch0/output.log")"
if [ "$rc1" -ne 0 ] && [ "$n1" -ge 1 ] && [ "$n0" -eq 0 ] && grep -q "main: done" "$work/switch0/output.log"; then
  echo "PASS"; exit 0
fi
echo "FAIL"; exit 1
