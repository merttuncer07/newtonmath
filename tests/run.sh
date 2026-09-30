#!/bin/sh
# Known-answer cases (tests/cases.txt), continuation of a root, and square roots checked against bc.
cd "$(dirname "$0")/.." || exit 2
NM=./newtonmath
fail=0; n=0
while IFS= read -r line; do
  case "$line" in ''|'#'*) continue;; esac
  stmt=${line%% ||| *}; want=${line#* ||| }
  got=$($NM -e "$stmt")
  n=$((n+1))
  if [ "$got" != "$want" ]; then echo "FAIL: $stmt"; echo "  want: $want"; echo "  got:  $got"; fail=$((fail+1)); fi
done < tests/cases.txt

# continuation: the second request continues from the stored root
out=$(printf 'let r = root of y^3 - 2y - 5 = 0 near 2\nr to 60 places\n' | $NM | tail -1)
n=$((n+1))
case "$out" in 2.094551481542326591482386540579302963857306105628239180304129\ *) ;; *) echo "FAIL: continuation: $out"; fail=$((fail+1));; esac

# a second route: bc computes the same square roots independently
if command -v bc >/dev/null 2>&1; then
  for k in 7 13 2 99991 123456789; do
    for places in 10 100 400; do
      mine=$($NM -e "sqrt($k) to $places places" | cut -d' ' -f1 | tr -d '.')
      theirs=$(printf 'scale=%d\nv=sqrt(%d)*10^%d+0.5\nscale=0\nv/1\n' $((places+10)) $k $places | bc | tr -d '\\\n')
      n=$((n+1))
      [ "$mine" = "$theirs" ] || { echo "FAIL: sqrt($k) to $places places differs from bc"; fail=$((fail+1)); }
    done
  done
fi
# values of series at numbers, against bc -l (a second, independent route)
if command -v bc >/dev/null 2>&1; then
  while IFS='|' read -r mine_expr bc_expr; do
    for places in 30 200; do
      mine=$(printf 'use prelude\n%s to %d places\n' "$mine_expr" $places | ./newtonmath | tail -1)
      case "$mine" in *"$places places guaranteed"*) ;; *) echo "FAIL: $mine_expr: $mine"; fail=$((fail+1)); continue;; esac
      mine=$(echo "$mine" | cut -d' ' -f1 | tr -d '.' | sed 's/^0*//')
      theirs=$(printf 'scale=%d\nv=(%s)*10^%d+0.5\nscale=0\nv/1\n' $((places+10)) "$bc_expr" $places | BC_LINE_LENGTH=0 bc -l)
      n=$((n+1))
      [ "$mine" = "$theirs" ] || { echo "FAIL: $mine_expr to $places places differs from bc"; fail=$((fail+1)); }
    done
  done <<'CASES'
exp(1)|e(1)
2 log1p(1/5) - log1p(-1/5) - log1p(-1/10)|l(2)
6 log1p(1/5) - 4 log1p(-1/5) - 3 log1p(-1/10)|l(10)
sin(1/2)|s(1/2)
cos(1/3)|c(1/3)
atan(1/5)|a(1/5)
6 asin(1/2)|4*a(1)
16 atan(1/5) - 4 atan(1/239)|4*a(1)
exp(-3/2)|e(-3/2)
CASES
fi

# series: the whole file against its reviewed output (every line was checked against Newton's texts or by hand)
n=$((n+1))
if ! ./newtonmath tests/series.nm | diff -u tests/series.out - ; then echo "FAIL: tests/series.nm"; fail=$((fail+1)); fi

echo "run.sh: $n checks, $fail failed"
[ $fail -eq 0 ]
