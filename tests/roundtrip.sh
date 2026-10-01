#!/bin/sh
# Slice 10: every printed exact quantity, rational function or matrix reads back as the same value.
# For each statement of the fixtures, the printed result R is fed back as (R) - (statement); it must be 0.
# Exempt: roots of letters (sqrt(x + 1) reads back as a series), series, verdict texts (isprime, factor), places, logs and arcs (log|..|, atan), roots r_k, verdict texts, solve and eliminate.
cd "$(dirname "$0")/.." || exit 2
NM=./newtonmath
S=7777777
total=0; bad=0
for f in tests/series.nm tests/newton.nm tests/irrational.nm tests/rules.nm tests/integration.nm tests/rational.nm tests/integral.nm tests/rootint.nm tests/euler.nm; do
  # pass 1: a sentinel after each statement splits the output by statement
  awk -v s=$S '!/^[ \t]*(#|$)/ { print; print s }' "$f" > /tmp/nm_rt1.$$
  $NM /tmp/nm_rt1.$$ > /tmp/nm_rt1o.$$ 2>&1
  # pass 2: the same statements, each candidate followed by its check, marked by a second sentinel
  awk -v s=$S -v k=0 -v i=0 'NR==FNR { if ($0 == (s "  [exact]")) { k++; next } if (seen[k]++) out[k] = out[k] "\n" $0; else out[k] = $0; next }
    /^[ \t]*(#|$)/ { next }
    { stmt = $0; o = out[i++]; print stmt
      if (o ~ /\n/ || o !~ /  \[exact/ || o ~ /log|atan|asin|sqrt\([^)]*[a-z]|O\(|r_[0-9]|where|any combination|error/) next
      if (stmt ~ /isprime|factor\(/ || stmt ~ /(^| )(to|solve|eliminate|use)( |$)/ || stmt ~ /^let [A-Za-z_]+[(\[]/) next
      body = o; sub(/  \[exact.*$/, "", body)
      if (stmt ~ /^let /) { name = stmt; sub(/^let /, "", name); sub(/ *=.*$/, "", name); sub(/^[^=]*= /, "", body); rhs = name }
      else rhs = "(" stmt ")"
      print "8888888"; print "(" body ") - " rhs }' /tmp/nm_rt1o.$$ "$f" > /tmp/nm_rt2.$$
  $NM /tmp/nm_rt2.$$ 2>&1 | awk -v f="$f" 'p { n++; if ($0 !~ /^(0|\[\[0(, 0)*\](, \[0(, 0)*\])*\])  \[exact/) { print "FAIL: round trip in " f ": " $0; b++ } p = 0; next }
    $0 == "8888888  [exact]" { p = 1 } END { print n + 0, b + 0 }' > /tmp/nm_rt3.$$
  grep '^FAIL' /tmp/nm_rt3.$$
  set -- $(tail -1 /tmp/nm_rt3.$$); total=$((total + $1)); bad=$((bad + $2))
done
rm -f /tmp/nm_rt*.$$
echo "roundtrip: $total checks, $bad failed"
[ $bad -eq 0 ]
