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

# Term order is independent of which letter was registered first.
for first in a x; do
  out=$(printf '%s\n(a+x)^2\n' "$first" | $NM | tail -1)
  n=$((n+1))
  [ "$out" = 'x^2 + 2ax + a^2  [exact]' ] || { echo "FAIL: term order after $first: $out"; fail=$((fail+1)); }
done
# Letters are single (aa is a*a; Slice 10), flowing letters lead; exact exponents even at total degree zero.
out=$(printf 'az\naa+az\n' | $NM | tail -1)
n=$((n+1))
[ "$out" = 'az + a^2  [exact]' ] || { echo "FAIL: multi-character term order: $out"; fail=$((fail+1)); }
out=$($NM -e 'a^(-1/2)*x^(1/2)+a^(1/2)*x^(-1/2)')
n=$((n+1))
[ "$out" = 'x^(1/2)/(a^(1/2)) + a^(1/2)/(x^(1/2))  [exact]' ] || { echo "FAIL: fractional term order: $out"; fail=$((fail+1)); }

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
exp(sqrt(2))|e(sqrt(2))
sin(sqrt(3)/2)|s(sqrt(3)/2)
log1p(sqrt(2) - 1)|l(sqrt(2))
atan(sqrt(3) - 1)|a(sqrt(3)-1)
sqrt(1 + sqrt(2))|sqrt(1+sqrt(2))
(3 + sqrt(5))/(1 + sqrt(2))|(3+sqrt(5))/(1+sqrt(2))
sin'(1/2)|c(1/2)
-cos''(1/3)|c(1/3)
CASES
  # values at complex points: real and imaginary parts, each against bc
  while IFS='|' read -r mine_expr bc_re bc_im; do
    places=60
    mine=$(printf 'use prelude\n%s to %d places\n' "$mine_expr" $places | ./newtonmath | tail -1)
    case "$mine" in *"$places places guaranteed"*) ;; *) echo "FAIL: $mine_expr: $mine"; fail=$((fail+1)); continue;; esac
    re=$(echo "$mine" | cut -d' ' -f1 | tr -d '.' | sed 's/^0*//')
    im=$(echo "$mine" | cut -d' ' -f3 | tr -d '.i' | sed 's/^0*//')
    sg=$(echo "$mine" | cut -d' ' -f2); [ "$sg" = "-" ] && im="-$im"
    for part in re im; do
      if [ $part = re ]; then e=$bc_re; got=$re; else e=$bc_im; got=$im; fi
      theirs=$(printf 'scale=%d\nv=(%s)*10^%d\nif (v<0) v=v-0.5 else v=v+0.5\nscale=0\nv/1\n' $((places+10)) "$e" $places | BC_LINE_LENGTH=0 bc -l)
      n=$((n+1))
      [ "$got" = "$theirs" ] || { echo "FAIL: $mine_expr ($part) differs from bc: $got vs $theirs"; fail=$((fail+1)); }
    done
  done <<'CASES'
exp(1 + i)|e(1)*c(1)|e(1)*s(1)
sin(1/2 + i/3)|s(1/2)*(e(1/3)+e(-1/3))/2|c(1/2)*(e(1/3)-e(-1/3))/2
log1p(i/2)|l(5/4)/2|a(1/2)
(1 + 2i) * exp(i)|c(1)-2*s(1)|s(1)+2*c(1)
CASES
fi

# Closed-form rational integrals: independent conic values at 60 places.
# Count a failed guarantee as a failed check too; round negatives symmetrically.
if command -v bc >/dev/null 2>&1; then
  while IFS='|' read -r mine_expr bc_expr; do
    places=60
    n=$((n+1))
    mine=$($NM -e "$mine_expr to $places places")
    case "$mine" in *"$places places guaranteed"*) ;; *) echo "FAIL: $mine_expr: $mine"; fail=$((fail+1)); continue;; esac
    mine=$(printf '%s\n' "$mine" | cut -d' ' -f1 | tr -d '.' | sed 's/^0*//; s/^-0*/-/')
    theirs=$(printf 'scale=%d\nv=(%s)*10^%d\nif (v<0) v=v-0.5 else v=v+0.5\nscale=0\nv/1\n' $((places+15)) "$bc_expr" $places | BC_LINE_LENGTH=0 bc -l | tr -d '\\\n')
    [ "$mine" = "$theirs" ] || { echo "FAIL: $mine_expr to $places places differs from bc: $mine vs $theirs"; fail=$((fail+1)); }
  done <<'CASES'
integral(1/(1+x),x,0,1)|l(2)
integral(1/(1+x^2),x,0,1)|a(1)
integral(1/(x^2+x+1),x,0,1)|4*a(1)*sqrt(3)/9
integral(x/(x^2+1),x,0,2)|l(5)/2
integral(1/(x^3+1),x,0,1)|l(2)/3+4*a(1)*sqrt(3)/9
integral(1/(1+x^2),x,0,2)|a(2)
integral(1/(1+x^2),x,0,100)|a(100)
integral(1/(1+x^2),x,0,-7)|-a(7)
integral(1/(1+x^2),x,-3/4,5/4)|a(5/4)+a(3/4)
integral(1/(1+x^2),x,1,0)|-a(1)
integral(1/x,x,1,1024)|l(1024)
integral(1/x,x,1,1/1024)|-l(1024)
integral(1/x,x,-4,-2)|-l(2)
integral(1/(1+x),x,1,0)|-l(2)
integral(1/(1+x^2),x,0,sqrt(3))|4*a(1)/3
integral(1/x,x,1,sqrt(2))|l(2)/2
integral(1/(x^2-2),x,2,3)|sqrt(2)/4*l((3-sqrt(2))*(2+sqrt(2))/((3+sqrt(2))*(2-sqrt(2))))
integral(1/(x^2+1)^2,x,0,1)|1/4+a(1)/2
integral((2x+3)/(x^2+2x+5),x,0,2)|l(13/5)+(a(3/2)-a(1/2))/2
integral(1/((x^2+1)*(x^2+2)),x,0,1)|a(1)-a(1/sqrt(2))/sqrt(2)
CASES
  # A user recipe can create the same printed name as a conic constant.
  # Its callback must not replace the private, mathematically defined atan(2).
  n=$((n+1))
  mine=$(printf "let atan = root of y' = 1, y(0) = 0\nlet T = atan(2 + x) to x^1\nintegral(1/(1+x^2),x,0,2) to 60 places\n" | $NM | tail -1)
  case "$mine" in
    *"60 places guaranteed"*)
      mine=$(printf '%s\n' "$mine" | cut -d' ' -f1 | tr -d '.')
      theirs=$(printf 'scale=75\nv=a(2)*10^60+0.5\nscale=0\nv/1\n' | BC_LINE_LENGTH=0 bc -l | tr -d '\\\n')
      [ "$mine" = "$theirs" ] || { echo "FAIL: conic atan(2) reused a user callback: $mine vs $theirs"; fail=$((fail+1)); };;
    *) echo "FAIL: conic atan(2) callback isolation: $mine"; fail=$((fail+1));;
  esac
else
  echo "SKIP: rational integral comparisons require bc -l"
fi

# The finite acceptance key, exact inverse checks, persistence and domain refusals.
n=$((n+1))
if ! ./newtonmath tests/integral.nm | diff -u tests/integral.out - ; then echo "FAIL: tests/integral.nm"; fail=$((fail+1)); fi

# 1000 places from the term rule, against bc (this took 17 s before the bounds were kept short; now under 1 s)
if command -v bc >/dev/null 2>&1; then
  for c in 'exp(1/3)|e(1/3)' 'atan(1/5)|a(1/5)'; do
    m=${c%%|*}; b=${c#*|}
    mine=$(printf 'use prelude\n%s to 1000 places\n' "$m" | ./newtonmath | tail -1)
    case "$mine" in *"1000 places guaranteed"*) ;; *) echo "FAIL: $m: $mine"; fail=$((fail+1)); continue;; esac
    mine=$(echo "$mine" | cut -d' ' -f1 | tr -d '.' | sed 's/^0*//')
    theirs=$(printf 'scale=1015\nv=(%s)*10^1000+0.5\nscale=0\nv/1\n' "$b" | BC_LINE_LENGTH=0 bc -l)
    n=$((n+1)); [ "$mine" = "$theirs" ] || { echo "FAIL: $m to 1000 places differs from bc"; fail=$((fail+1)); }
  done
fi

# series: the whole file against its reviewed output (every line was checked against Newton's texts or by hand)
n=$((n+1))
if ! ./newtonmath tests/series.nm | diff -u tests/series.out - ; then echo "FAIL: tests/series.nm"; fail=$((fail+1)); fi

# letters and the parallelogram: Newton's own examples (see the comments in the file)
n=$((n+1))
if ! ./newtonmath tests/newton.nm | diff -u tests/newton.out - ; then echo "FAIL: tests/newton.nm"; fail=$((fail+1)); fi

# surds and i
n=$((n+1))
if ! ./newtonmath tests/irrational.nm | diff -u tests/irrational.out - ; then echo "FAIL: tests/irrational.nm"; fail=$((fail+1)); fi

# the four modules reached from the language: systems, matrices, whole numbers, complex points, moved series
n=$((n+1))
if ! ./newtonmath tests/integration.nm | diff -u tests/integration.out - ; then echo "FAIL: tests/integration.nm"; fail=$((fail+1)); fi

# Rational functions and matrices: reviewed answer key, persistence, domains and refusals.
n=$((n+1))
if ! ./newtonmath tests/rational.nm | diff -u tests/rational.out - ; then echo "FAIL: tests/rational.nm"; fail=$((fail+1)); fi
# Fresh processes see the letters in opposite orders; both expressions have one normal form.
for first in a x; do
  for expr in '(x-a)/(a^2-x^2)' '-1/(x+a)'; do
    out=$(printf '%s\n%s\n' "$first" "$expr" | $NM | tail -1)
    n=$((n+1))
    [ "$out" = '-1/(x + a)  [exact]' ] || { echo "FAIL: quotient order after $first: $out"; fail=$((fail+1)); }
  done
done
for first in a d; do
  out=$(printf '%s\nrank([[a,b],[c,d]])\n' "$first" | $NM | tail -1)
  n=$((n+1))
  [ "$out" = '2  [exact, for general a, b, c, d; rank drops where ad - bc = 0]' ] || { echo "FAIL: rank order after $first: $out"; fail=$((fail+1)); }
done

# Slice 11: areas under roots, the whole file, and definite values against bc -l to 60 places
n=$((n+1))
if ! ./newtonmath tests/rootint.nm | diff -u tests/rootint.out - ; then echo "FAIL: tests/rootint.nm"; fail=$((fail+1)); fi
n=$((n+1))
if ! ./newtonmath tests/euler.nm | diff -u tests/euler.out - ; then echo "FAIL: tests/euler.nm"; fail=$((fail+1)); fi
if command -v bc >/dev/null 2>&1; then
  while IFS='|' read -r mine_expr bc_expr; do
    places=60
    mine=$($NM -e "$mine_expr to $places places")
    case "$mine" in *"$places places guaranteed"*) ;; *) echo "FAIL: $mine_expr: $mine"; fail=$((fail+1)); continue;; esac
    mine=$(echo "$mine" | cut -d' ' -f1 | tr -d '.' | sed 's/^0*//')
    theirs=$(printf 'scale=%d\nv=(%s)*10^%d+0.5\nscale=0\nv/1\n' $((places+15)) "$bc_expr" $places | BC_LINE_LENGTH=0 bc -l)
    n=$((n+1))
    [ "$mine" = "$theirs" ] || { echo "FAIL: $mine_expr differs from bc"; fail=$((fail+1)); }
  done <<'CASES'
integral(sqrt(1 - x^2), x, 0, 1)|a(1)
integral(1/sqrt(1 - x^2), x, 0, 1/2)|4*a(1)/6
integral(1/sqrt(x^2 + 1), x, 0, 1)|l(1+sqrt(2))
integral(sqrt(x^2 - 1), x, 1, 2)|sqrt(3)-l(2+sqrt(3))/2
integral(x^2 sqrt(x^2 - 1), x, 1, 2)|(16*sqrt(3)-2*sqrt(3)-l(2+sqrt(3)))/8
integral(1/(x sqrt(x + 1)), x, 1, 3)|l((2-1)/(2+1))-l((sqrt(2)-1)/(sqrt(2)+1))
integral((x + 1)/sqrt(x^2 + 2x + 5), x, 0, 1)|sqrt(8)-sqrt(5)
integral(x^3/sqrt(1 - x^2), x, 0, 1/2)|2/3-(1/4+2)*sqrt(3/4)/3
integral(1/sqrt(2x - x^2), x, 1/2, 1)|4*a(1)/6
integral(1/(x sqrt(x^2 + 1)), x, 1, 2)|l(1+sqrt(2))-l((1+sqrt(5))/2)
integral(1/(x^2 sqrt(x^2 + 1)), x, 1, 2)|sqrt(2)-sqrt(5)/2
integral(1/(x sqrt(x^2 - 2)), x, 2, 3)|(a(sqrt(7)/sqrt(2))-a(1))/sqrt(2)
CASES
fi

# Slice 14: factors with no split over Q; definite values against independent quadrature (sympy, 25 places)
n=$((n+1))
if ! ./newtonmath tests/realsplit.nm | diff -u tests/realsplit.out - ; then echo "FAIL: tests/realsplit.nm"; fail=$((fail+1)); fi
while IFS='|' read -r f a b want; do
  got=$($NM -e "integral($f, x, $a, $b) to 25 places" | cut -d' ' -f1)
  n=$((n+1)); [ "$got" = "$want" ] || { echo "FAIL: integral($f, x, $a, $b): $got, quadrature $want"; fail=$((fail+1)); }
done <<'CASES'
1/(x^3 + x + 1)|0|1|0.6303193224124080140667863
1/(x^3 - 3x + 1)|-1|0|0.4895517625859923877210842
x/(x^3 - 2)|0|1|-0.3251555784294800194916487
1/(x^4 + 1)|0|1|0.8669729873399110375739952
1/(x^4 + 2)|0|2|0.6208412393004786044234699
x^2/(x^4 + x + 1)|0|1|0.1616404362034916102203722
(x^2 + 1)/(x^4 - 10x^2 + 1)|1|2|-0.2088101520505993212096815
1/(x^3 + 2x + 5)|0|3|0.2923682218213270459150485
CASES

# sparse letters: one long session (Slices 11 and 12 together, and 240 surds and named constants) stays exact
n=$((n+1))
if ! cat tests/rootint.nm tests/euler.nm | ./newtonmath | diff -q - tests/rootint.out tests/euler.out >/dev/null 2>&1; then
  cat tests/rootint.out tests/euler.out > /tmp/nm_both.$$
  cat tests/rootint.nm tests/euler.nm | ./newtonmath | diff -q - /tmp/nm_both.$$ >/dev/null || { echo "FAIL: one long session"; fail=$((fail+1)); }
  rm -f /tmp/nm_both.$$
fi
out=$(k=1; while [ $k -le 120 ]; do echo "sqrt($((2*k+1))) + log($((k+1)))"; k=$((k+1)); done; echo '(a+b+c+d+e+f+g+h+j+k+l+m)^2 - (m+l+k+j+h+g+f+e+d+c+b+a)^2')
out=$(echo "$out" | ./newtonmath | tail -1)
n=$((n+1)); [ "$out" = '0  [exact]' ] || { echo "FAIL: many letters: $out"; fail=$((fail+1)); }

# rules, sequences, sums and cases
n=$((n+1))
if ! ./newtonmath tests/rules.nm | diff -u tests/rules.out - ; then echo "FAIL: tests/rules.nm"; fail=$((fail+1)); fi
# a deep rule and a long table, against bc
if command -v bc >/dev/null 2>&1; then
  mine=$(printf 'let fact(n) = 1 if n = 0, n * fact(n - 1) otherwise\nfact(300)\n' | ./newtonmath | tail -1 | cut -d' ' -f1)
  theirs=$(printf 'p=1\nfor(i=1;i<=300;i++) p*=i\np\n' | BC_LINE_LENGTH=0 bc)
  n=$((n+1)); [ "$mine" = "$theirs" ] || { echo "FAIL: fact(300) differs from bc"; fail=$((fail+1)); }
  mine=$(printf 'let F[0] = 0, F[1] = 1, F[k] = F[k - 1] + F[k - 2]\nF[2000]\n' | ./newtonmath | tail -1 | cut -d' ' -f1)
  theirs=$(printf 'a=0\nb=1\nfor(i=0;i<2000;i++){t=a+b;a=b;b=t}\na\n' | BC_LINE_LENGTH=0 bc)
  n=$((n+1)); [ "$mine" = "$theirs" ] || { echo "FAIL: F[2000] differs from bc"; fail=$((fail+1)); }
fi

# show: Newton's own table for y^3 - 2y - 5 (De Analysi, NATP00204), and the same answer as without show.
out=$($NM -e 'show root of y^3 - 2y - 5 = 0 near 2 to 30 places')
n=$((n+1))
echo "$out" | grep -qx '    y = 2 + p:   p^3 + 6p^2 + 10p - 1 = 0,   p = 0.1' || { echo "FAIL: show, first row"; fail=$((fail+1)); }
n=$((n+1))
echo "$out" | grep -qx '    y = 2.1 + p:   p^3 + 6.3p^2 + 11.23p + 0.061 = 0,   p = -0.0054' || { echo "FAIL: show, second row"; fail=$((fail+1)); }
n=$((n+1))
[ "$(echo "$out" | tail -1)" = "$($NM -e 'root of y^3 - 2y - 5 = 0 near 2 to 30 places')" ] || { echo "FAIL: show changes the answer"; fail=$((fail+1)); }
out=$($NM -e 'show 2 + 3' | head -1)
n=$((n+1))
[ "$out" = '  (no steps are written for this kind of statement yet)' ] || { echo "FAIL: show without steps: $out"; fail=$((fail+1)); }

# The libraries are compiled into the program: it runs alone, from any directory.
tmp=$(mktemp -d); cp $NM "$tmp/"
out=$(cd "$tmp" && printf 'use prelude\nexp(1) to 30 places\n' | ./newtonmath | tail -1)
n=$((n+1))
[ "$out" = '2.718281828459045235360287471353  [bounded: 30 places guaranteed]' ] || { echo "FAIL: built-in prelude: $out"; fail=$((fail+1)); }
echo 'let w = 7' > "$tmp/mine.nm"
out=$(cd "$tmp" && printf 'use mine\nw\n' | ./newtonmath | tail -1)
n=$((n+1))
[ "$out" = '7  [exact]' ] || { echo "FAIL: library from this directory: $out"; fail=$((fail+1)); }
rm -rf "$tmp"

echo "run.sh: $n checks, $fail failed"
[ $fail -eq 0 ]
