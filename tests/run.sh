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
echo "run.sh: $n checks, $fail failed"
[ $fail -eq 0 ]
