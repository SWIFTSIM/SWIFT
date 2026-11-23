#!/bin/bash
#
# Key arrival test of the GeajayCooling parameters.
#
# Needs a SWIFT built with --with-cooling=geajay_3, the program testGeajayKeys
# (make -C tests testGeajayKeys) and four Cloudy tables. Not part of "make
# check". Usage:
#
#   tests/testGeajayKeys.sh <testGeajayKeys> <fixture dir> <work dir> [cases file]
#
# The variable BASE_EXTRA adds lines (separated by ";") to the baseline file of
# every case, for example BASE_EXTRA="background_mode: per_band" to send the
# options of the cases through the per_band entry point of the library.
#
# The fixture dir holds CloudyData_UVB=HM2012.h5 (the table of the example
# CoolingBox), colA.h5 (a Cloudy table with the column shielding), g0file.h5 (a
# Cloudy table with an interstellar radiation field) and hm12thin.h5 (an
# optically thin UV background table). Their content is read, so they must be
# valid tables.
#
# Each line of the cases file is one test:
#   id|key|base lines|test lines|expected|hint
# The program reads a baseline file and a test file that differs by the test
# lines. The test passes when every expected line is in the output of the test
# file and not in the output of the baseline file. An expected "~MOD name" line
# only needs a different value. An expected "!text" makes the test file a
# refusal: the program must stop with text in its output. An expected "=line"
# only needs the line in the output of the test file (a default that must not
# move). An expected "@key" is a refusal that suggests a value ("Set
# GeajayCooling:key to <value>,"): the test file must stop, and a third file
# with the suggested value must start. The lines of a file are separated by ";", a later line replaces an
# earlier key and "-key" drops one. Exit code 0: every case passes.

here=$(dirname "$(realpath "$0")")
bin=$(realpath "$1")
fixture=$(realpath "$2")
mkdir -p "$3" || exit 2
out=$(realpath "$3")
cases=$(realpath "${4:-$here/testGeajayKeys.cases}")
head=$here/testGeajayKeys.yml
unset LD_LIBRARY_PATH

cd "$out" || exit 2
for f in CloudyData_UVB=HM2012.h5 colA.h5 g0file.h5 hm12thin.h5; do
  ln -sf "$fixture/$f" "$f"
done
ln -sf "$fixture/CloudyData_UVB=HM2012.h5" alt_table.h5
ln -sf "$fixture/CloudyData_UVB=HM2012.h5" uvrates.h5
ln -sf "$fixture/colA.h5" colB.h5
ln -sf "$fixture/g0file.h5" g0file2.h5

# mk <file> <lines>
mk() {
  local IFS=';' l spec=""
  for l in $2; do
    [ -z "$l" ] && continue
    spec="$spec$l"$'\n'
  done
  awk -v spec="$spec" '
    function key(s) { sub(/^ +/, "", s); sub(/:.*/, "", s); return s }
    BEGIN { n = split(spec, a, "\n")
      for (i = 1; i <= n; i++) { l = a[i]; if (l == "") continue
        if (substr(l, 1, 1) == "-") { delete over[substr(l, 2)]; gone[substr(l, 2)] = 1 }
        else { k = key(l); over[k] = "  " l; delete gone[k]; if (!(k in seen)) { seen[k] = 1; ord[++m] = k } } } }
    /^  [a-zA-Z]/ && inblock { k = key($0); if (k in gone) next; if (k in over) { print over[k]; used[k] = 1; next } }
    /^[A-Za-z]/ { if (inblock) flush(); inblock = ($0 == "GeajayCooling:") }
    { print }
    function flush(  i) { for (i = 1; i <= m; i++) if ((ord[i] in over) && !(ord[i] in used)) print over[ord[i]] }
    END { if (inblock) flush() }' "$head" > "$1"
}

fails=0
n=0
while IFS='|' read -r id key base test expected hint; do
  [[ "$id" == \#* || -z "$id" ]] && continue
  n=$((n + 1))
  d="$out/$id"
  mkdir -p "$d"
  base="$base${BASE_EXTRA:+;$BASE_EXTRA}"
  mk "$d/base.yml" "$base"
  mk "$d/test.yml" "$base;$test"
  "$bin" "$d/base.yml" > "$d/base.out" 2>&1
  sb=$?
  "$bin" "$d/test.yml" > "$d/test.out" 2>&1
  st=$?
  verdict=PASS
  evid=""
  if [ $sb -ne 0 ] || ! grep -q KEYPASS_DONE "$d/base.out"; then
    verdict="FAIL(baseline exit $sb)"
  fi
  if [[ "$expected" == "!"* ]]; then
    text=${expected#!}
    if [ $st -ne 0 ] && grep -qF -- "$text" "$d/test.out" &&
      ! grep -qF -- "$text" "$d/base.out"; then
      evid=" [refused: $text]"
    else
      verdict="FAIL(no refusal with: $text; exit $st)"
    fi
  elif [[ "$expected" == "@"* ]]; then
    name=${expected#@}
    value=$(sed -n "s/.*GeajayCooling:$name to \([^,]*\),.*/\1/p" "$d/test.out" | head -n 1)
    if [ $st -eq 0 ] || [ -z "$value" ]; then
      verdict="FAIL(no suggested value for $name; exit $st)"
    else
      mk "$d/fed.yml" "$base;$test;$name: $value"
      "$bin" "$d/fed.yml" > "$d/fed.out" 2>&1
      sf=$?
      if [ $sf -eq 0 ] && grep -q KEYPASS_DONE "$d/fed.out"; then
        evid=" [suggested $name: $value starts]"
      else
        verdict="FAIL(suggested $name: $value stops; exit $sf: $(tail -n 2 "$d/fed.out" | tr '\n' ' ' | cut -c1-200))"
      fi
    fi
  elif [ $st -ne 0 ] || ! grep -q KEYPASS_DONE "$d/test.out"; then
    verdict="FAIL(test exit $st: $(tail -n 2 "$d/test.out" | tr '\n' ' ' | cut -c1-300))"
  else
    IFS=$'\n' read -r -d '' -a exp < <(
      printf '%s' "${expected//&&/$'\n'}"
      printf '\0'
    )
    for e in "${exp[@]}"; do
      if [[ "$e" == "~"* ]]; then
        name=${e#\~}
        a=$(grep -m1 "^$name=" "$d/base.out")
        b=$(grep -m1 "^$name=" "$d/test.out")
        if [ -n "$b" ] && [ "$a" != "$b" ]; then
          evid="$evid [$a -> $b]"
        else
          verdict="FAIL(no change in $name)"
        fi
      elif [[ "$e" == "="* ]]; then
        if grep -qxF -- "${e#=}" "$d/test.out"; then
          evid="$evid [${e#=}]"
        else
          verdict="FAIL(missing: ${e#=})"
        fi
      else
        if grep -qxF "$e" "$d/test.out" && ! grep -qxF "$e" "$d/base.out"; then
          evid="$evid [$e]"
        else
          verdict="FAIL(missing: $e)"
        fi
      fi
    done
  fi
  [[ "$verdict" == PASS ]] || fails=$((fails + 1))
  printf '%s\t%s\t%s\t%s\n' "$id" "$key" "$verdict" "$evid"
done < "$cases"
echo "cases: $n failures: $fails"
[ $fails -eq 0 ]
