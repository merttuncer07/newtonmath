# Slice 8 verification record

Verified on 2026-10-01, macOS arm64, Apple Clang 21, C11, normal `-O2` build. Implementation commits:
`514b63f` (isolated module) and `20090ef` (language/matrices/tests). No runtime dependencies were added.

## Acceptance answer key

The following are the actual outputs of the first 16 statements in tests/rational.nm. Ordinary finite results
carry `[exact]`. Matrix qualifications follow the table. Products in input use `*`; juxtaposed letters in the
output follow the existing printer. Equivalent sign choices follow the fixed leading-term order.

| Input | Printed result |
|---|---|
| `(x^2-a^2)/(x-a)` | `a + x` |
| `(x^3-a^3)/(x-a)` | `a^2 + ax + x^2` |
| `(x^4-a^4)/(x^2-a^2)` | `a^2 + x^2` |
| `(x^2-1)/(x^2+2x+1)` | `(x - 1)/(x + 1)` |
| `(x^2+2*a*x+a^2)/(x^2-a^2)` | `(-a - x)/(a - x)` |
| `1/(x+a)+1/(x-a)` | `-2x/(a^2 - x^2)` |
| `a/(a-b)+b/(b-a)` | `1` |
| `(a/b)/(c/d)` | `ad/bc` |
| `(a*x+a)/(b*x+b)` | `a/b` |
| `gcd(6*a*x+6a,4x+4)` | `x + 1` |
| `inverse([[a,b],[c,d]])` | `[[d/(ad - bc), -b/(ad - bc)], [-c/(ad - bc), a/(ad - bc)]]` |
| `inverse([[1,a],[a,1]])` | `[[-1/(a^2 - 1), a/(a^2 - 1)], [a/(a^2 - 1), -1/(a^2 - 1)]]` |
| `linsolve([[a,1],[1,a]],[1,0])` | `[[a/(a^2 - 1)], [-1/(a^2 - 1)]]` |
| `rank([[a,b],[c,d]])` | `2` |
| `a^2/(b+x)` | `a^2/(b + x)` |
| `a^2/(b+x) to x^3` | `a^2/b - a^2x/(b^2) + a^2x^2/(b^3) - a^2x^3/(b^4) + O(x^4)` |

Full matrix verdicts:

- General inverse: `[exact, for general a, b, c, d; displayed formula undefined where ad - bc = 0]`.
- Symmetric inverse: `[exact, for general a; displayed formula undefined where a^2 - 1 = 0]`.
- Linear solution: `[exact, for general a; rank drops where a^2 - 1 = 0; displayed formula undefined where a^2 - 1 = 0]`.
- Rank: `[exact, for general a, b, c, d; rank drops where ad - bc = 0]`.

Checks against the supplied independent key: difference-of-powers factorizations give the first three rows;
`(x+1)^2` gives row 4; `(x+a)^2/(x^2-a^2)` cancels x+a and agrees with row 5 after negating numerator and
denominator. Row 6 has common denominator x^2-a^2 and numerator 2x; the monic convention negates both. Row 7
is `(a-b)/(a-b)=1`; rows 8-9 cancel the inner denominators and the common x+1. Over Q, the content factor a
of the first gcd input is absent from the second input, leaving x+1. The 2x2 adjugate identities verify the
matrix rows. Newton's geometric division of aa/(b+x) supplies the last row.

## Check counts

`make test` reports **10,564 checks**, all passing. Fixture comparisons count as one shell check each; they
are not individually added again to this total. The Makefile runs each unit executable twice, once silently
for failure handling and once for the summary; the following counts describe one pass per executable.

| Suite | Reported checks | Result |
|---|---:|---|
| tests/run.sh | 104 | 0 failed |
| tests/check_z | 9065 | 0 failed |
| unit arith | 19 | 0 failed |
| unit cplx | 13 | 0 failed |
| unit elim | 9 | 0 failed |
| unit linalg | 13 | 0 failed |
| unit ratfun | 788 | 0 failed |
| unit ratmat | 553 | 0 failed |

Breakdown of the new module tests:

| Kind | ratfun | ratmat |
|---|---:|---:|
| Symbolic multiplication back | 11 | 10 |
| Lowest-terms gcd checks | 30 | 25 |
| Positive leading denominator | 11 | included in normal-form cases |
| Exact printed normal forms | 19 | covered in language fixture |
| Gcd/content/exact-division edge cases | 8 | — |
| Refusal with reason | 4 | 1 |
| Existing surd inverse compatibility | 1 | covered by old linalg tests |
| Exact rank/domain/basis-pole conditions | — | 5 |
| Independent Q comparisons | 704 | 512 |

- Ratfun: 64 accepted sets of small rational values, 24 singular candidates skipped; fixed seed 104729. All
  11 original expressions are evaluated directly with Q, independently of R/poly_gcd, at each accepted set.
- Ratmat: 64 accepted sets, 3 singular candidates skipped; fixed seed 65537. Four inverse entries are compared
  with the Q adjugate formula and four entries of the product with the original numeric matrix are checked.
  The first case has a=0 with nonzero determinant, so a pivot choice cannot masquerade as a singular locus.
- Ratfun was also run in a fresh process registering x before a: the same 788 checks passed. The language
  runner adds four quotient normalization checks under both registration orders and two rank-order checks.
- tests/rational.nm has **67 statements**, including **six expected refusals**: fractional quotient powers,
  mixed letter/surd denominator, division by zero, identically singular inverse, and two unsupported prime
  syntax uses on a bound rational quantity. The explicit surd-series case after the refusals still succeeds.
- Production gcd calls divide both original inputs exactly and check division by multiplication. Each
  nontrivial quotient normalization cross-multiplies with the original numerator and denominator. These
  internal checks execute on every operation and are not separately counted as regression assertions.

The fixed-seed numeric checks are regression evidence, not randomized proofs. Exactness comes from the
rational arithmetic and deterministic polynomial algorithms; no new result is labelled `[proved]`.

## Measured times and compiler coverage

| Command / suite | Time |
|---|---:|
| `make test`, after final code edit (includes rebuilding all test executables) | 12.95 s wall |
| `tests/unit/ratfun_test` | 9.14 ms CPU |
| `tests/unit/ratfun_test x-first` | 8.24 ms CPU |
| `tests/unit/ratmat_test` | 134.23 ms CPU |
| 67-statement rational language fixture, normal binary | 36.69 ms wall |
| Same fixture under AddressSanitizer + UndefinedBehaviorSanitizer | 481.52 ms wall |

The sanitizer run exactly matched the reviewed fixture and produced no diagnostics. Leak detection was
explicitly disabled because the language's permanent arena lives until process exit; invalid memory accesses
and undefined behavior remained enabled with halt-on-error.

Both `/usr/bin/clang` and `/usr/bin/gcc` identify as **Apple Clang 21.0.0**. Building through the gcc driver
produced identical output for the entire new fixture. **An independent GNU GCC compiler was not available**;
this is driver equivalence, not two-compiler coverage. No compiler was installed for this task.

## Every changed pre-existing .out line

Only **tests/integration.out:70** changes (input tests/integration.nm:65):

Before:

```text
g = 1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + O(x^9)
```

After:

```text
g = -1/(x - 1)  [exact]
```

Reason: the definition is `let g = 1/(1-x)`. The numerator and denominator of its new output multiply back
to that same function, while the old line was its expansion through degree 8. D1 requires the finite answer
by default. The immediately following outputs remain `3/2 [exact]` at x=1/3 and the explicit moved series
`-1 + x - x^2 + x^3 + O(x^4)` at x=2+x. All other old fixture lines are unchanged.

## Deliberate scope limits

- Fractional exponents in finite quotients and mixed surd/ordinary-letter denominators are refused; polynomial
  gcd is over Q. Surd numerators over ordinary-letter denominators are also outside this v1 coefficient field.
- Series coefficients remain C; rational-function coefficients needing division by a sum such as a+b are not
  added to the series engine. Existing explicit surd expansions continue to work.
- Generic rank/solution conditions are stated, but exceptional parameter cases are not all solved or given
  alternate bases automatically. Maximal-minor enumeration can grow combinatorially on large matrices.
- Partial fractions, polynomial factorization and closed-form integration of rational functions remain next
  steps. Prime syntax on bound rational quantities is refused; use d/dx and then substitute.
- Cancellation is equality in the rational-function field, not bookkeeping for holes of the original expression.
