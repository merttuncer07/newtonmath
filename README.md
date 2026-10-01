# newtonmath

A language only for mathematical computation. It is written in C with no dependencies. Build with `make`;
test with `make test`.

    ./newtonmath                    # line by line
    ./newtonmath file.nm            # run a file
    ./newtonmath -e "sqrt(2) to 30 places"

Every result carries its verdict:
- `[exact]`: the value itself.
- `[certified: sign change of ..., k places]`: a root whose bracket was confirmed by substituting both ends.
- `[bounded: k places guaranteed]`: ball arithmetic; the printed places are guaranteed.
- `error: ...`: the statement could not be carried out, with the reason.

Statements:
- `let r = root of y^3 - 2y - 5 = 0 near 2`, then `r to 50 places`. The second request continues from the first.
- Plain ASCII is canonical. `√ × · − ² ³` are accepted as aliases.
- Expression operands are evaluated from left to right.

## Series (slice 2)

One engine, Newton's resolution, works for numbers and for series. The same `root of` resolves an algebraic
equation or a fluxional (differential) one from its start.

    sqrt(1 - x^2) to x^8                              1 - x^2/2 - x^4/8 - x^6/16 - 5x^8/128 + O(x^9)
    root of y^2 = 1 - x^2, y(0) = 1 to x^8            the same series, as the root of its equation
    root of y' = 1 - 3x + y + x^2 + x*y, y(0) = 0     Newton's Methodus, Problem 2
    root of y'' = -y, y(0) = 0, y'(0) = 1             sin
    integral(1/(1 + x)) to x^6                        log(1 + x)
    d/dx (x^3 - 2x)                                   3x^2 - 2  [exact]

Every series is substituted back into its equation, or raised back to its power, before it is shown.
Named functions are not built in. `use prelude` loads lib/prelude.nm, where exp, sin, cos, log1p, atan and asin
are defined by their equations, in the language itself. A defined series can be substituted: `exp(x^2)`.

## Values (slice 3)

A function defined by a linear fluxional equation with polynomial coefficients has a term rule: each coefficient
follows from the ones before. The language derives that rule from the equation, sums the terms, and bounds the
rest from the rule, so the places it prints are guaranteed. Arguments are exact numbers. As Newton did, you choose
small arguments and join them by exact relations:

    use prelude
    2 log1p(1/5) - log1p(-1/5) - log1p(-1/10) to 30 places     log 2 (1.2 * 1.2 / (0.8 * 0.9) = 2)
    16 atan(1/5) - 4 atan(1/239) to 50 places                  pi (Machin)
    exp(1) to 1000 places

Near the edge of convergence, for example log1p(1), the language says the terms shrink too slowly there.

## Letters and the parallelogram (slice 4)

Letters are given quantities, following Newton's convention: a, b, c are given and x, z flow. Exponents may be
fractions or negative (a^(1/2), x^(-2/5)). Quotients in ordinary letters with whole exponents use rational
functions (slice 8 below). `to x^8` explicitly requests a series and names its letter.

Finite sums print in descending total degree. Equal degrees are ordered by descending exponents of the
alphabetically ordered letter names, independently of when the letters first appeared:

    (a + x)^2                                  a^2 + 2ax + x^2  [exact]

    sqrt(a^2 + x^2) to x^6                     a + x^2/(2a) - x^4/(8a^3) + x^6/(16a^5) + O(x^7)
    root of y^3 + a^2 y + a x y - 2a^3 - x^3 = 0 for y
        2 starts; choose one with 'starting y = ...':  y = a  and two that are not rational numbers
    root of y^3 + a^2 y + a x y - 2a^3 - x^3 = 0 for y starting y = a to x^4
        a - x/4 + x^2/(64a) + 131x^3/(512a^2) + 509x^4/(16384a^3) + O(x^5)     (De analysi)

The parallelogram finds the first term of every root. A start that is a multiple root is refined by applying it
again. The finished root is substituted back into the equation.

## Rules, sequences and sums (slice 5)

Newton wrote repetition as a rule: each term from the one before, and tables built forward. The language has
rules, sequences defined by their earlier terms, definition by cases, and sums:

    let fact(n) = 1 if n = 0, n * fact(n - 1) otherwise
    let fib[0] = 0, fib[1] = 1, fib[k] = fib[k - 1] + fib[k - 2]
    let binom[0] = 1, binom[k] = (m - k + 1)/k * binom[k - 1]       Newton's binomial rule
    sum(binom[k] x^k for k = 0 to 8) - sqrt(1 + x) to x^8             O(x^9): agrees with the engine

A sequence keeps its table and is built forward, so fib[2000] needs no deep recursion. A definition by
cases compares exact numbers. Rules work on numbers, letters and series. This is the way toward writing the
language's own methods in the language (see docs/DECISIONS.md).

## Surds and i (slice 6)

Following Newton ("pro singulis pono totidem literas"), an irrational number is a letter with its own equation:
sqrt(2) satisfies r^2 = 2. Arithmetic stays exact; powers are reduced by the equation; surds are cleared from
denominators:

    (1 + sqrt(2))^2              2sqrt(2) + 3  [exact]
    1/(1 + 2^(1/3))              (2^(1/3))^2/3 - 2^(1/3)/3 + 1/3  [exact]
    let r = root of y^3 - 2y - 5 = 0 near 2
    1/r                          r^2/5 - 2/5  [exact]
    (3 + 4i)/(1 - 2i)            -1 + 2i  (i is the square root of -1)
    sqrt(2 + x) to x^3           sqrt(2) + sqrt(2)x/4 - sqrt(2)x^2/32 + sqrt(2)x^3/128 + O(x^4)
    exp(sqrt(2)) to 30 places    4.113250378782927517173581815140  (use prelude)

A named root whose equation factors is reduced to the factor its root satisfies. The parallelogram follows
starts with surds: starting y = sqrt(2) a^(1/2) x^(1/2).

## Systems, matrices, whole numbers, complex points (slice 7)

Four parts built apart (src/elim.c, linalg.c, cplx.c, arith.c, each with its own tests in tests/unit) and then
joined to the language in one pass:

    solve x^2 + y^2 = 5, x - y = 1 for x, y      x = -1, y = -2 / x = 2, y = 1  [each put back: exact]
    solve x^3 - x - 1 for x                      x = r_1, where r_1 = 1.3247...  [certified]
    eliminate y from x^2 + y^2 = 1, y = x^2      x^4 + x^2 - 1  [exact]
    let A = [[1, 2], [3, 4]]
    inverse(A), det(A), A^-2, transpose(A), rank(...), nullspace(...), linsolve(A, [3, 1])
    charpoly(A)                                  t^2 - 5t - 2
    eigenvalues(A)                               [[sqrt(33)/2 + 5/2, -sqrt(33)/2 + 5/2]]
    factor(2^64 + 1)                             274177 * 67280421310721  [multiplied back; every factor proved prime]
    isprime(2^61 - 1)                            true  [proved]
    gcd, lcm, mod, powmod, invmod, divisors, sigma, phi, nextprime
    exp(1 + i)                                   1.46869393991588515714 + 2.28735528717884239121i  (use prelude)
    sin'(1/2)                                    0.87758256189037271612  (a fluxion at a point)
    sin(1/2 + x) to x^3                          sin(1/2) + sin'(1/2)x - sin(1/2)x^2/2 - sin'(1/2)x^3/6 + O(x^4)

`[a, b]` is a column, `[[..], [..]]` a matrix. A probable prime is labelled probable, never proved. A series moved
to a point keeps the values at that point as letters (sin(1/2), sin'(1/2)), each known by its value.

Not yet:
- division by quantities with several different surds nested in each other's equations;
- arithmetic on fractional-power results;
- moving a series to a complex point.


## Rational functions (slice 8)

A finite quotient is exact by default. Ask for `to x^N` to expand it as a series. Ordinary letters with whole
exponents and rational coefficients form a field; recursive polynomial gcd reduces every fraction. The gcd
and denominator are monic in the fixed term order described above, so signs can differ from handwritten forms.
A unit denominator collapses to an ordinary number or polynomial.

    (x^2 - a^2)/(x - a)                    a + x  [exact]
    1/(x + a) + 1/(x - a)                  -2x/(a^2 - x^2)  [exact]
    a/(a - b) + b/(b - a)                  1  [exact]
    (a/b)/(c/d)                            ad/bc  [exact]
    gcd(6*a*x + 6*a, 4*x + 4)              x + 1  [exact]
    a^2/(b + x)                           a^2/(b + x)  [exact]
    a^2/(b + x) to x^3                    a^2/b - a^2x/(b^2) + a^2x^2/(b^3) - a^2x^3/(b^4) + O(x^4)

Write products of letters with `*` or spaces: `a*x` is a product, whereas `ax` is a name. Negative whole powers
are cleared into numerator and denominator. Polynomial gcd accepts nonnegative whole powers. Integer-only
`gcd(12,18)` retains the integer answer 6; a polynomial argument selects polynomial gcd over Q.

    let g = 1/(1 - x)
    g(1/3)                                3/2  [exact]
    g(2 + x) to x^3                       -1 + x - x^2 + x^3 + O(x^4)
    d/dx (1/(1 - x))                      1/(x^2 - 2x + 1)  [exact]

Rational functions can be stored, passed to rules, used in sequences and compared for exact equality. A stored
quotient in one ordinary letter accepts exact substitution; for multiple letters, define an explicit rule.
Equality is equality in the rational-function field: cancellation forgets removed holes, so `p/p` becomes 1.
The quotient operation still refuses an identically zero denominator.

Matrices accept rational-function entries. Inverse, linsolve, rank and nullspace work over this field:

    inverse([[a, b], [c, d]])
        [[d/(ad - bc), -b/(ad - bc)], [-c/(ad - bc), a/(ad - bc)]]
        [exact, for general a, b, c, d; displayed formula undefined where ad - bc = 0]
    linsolve([[a, 1], [1, a]], [1, 0])
        [[a/(a^2 - 1)], [-1/(a^2 - 1)]]
        [exact, for general a; rank drops where a^2 - 1 = 0; displayed formula undefined where a^2 - 1 = 0]
    rank([[a, b], [c, d]])
        2  [exact, for general a, b, c, d; rank drops where ad - bc = 0]
    nullspace([[a, b]])
        [[-b/a], [1]]
        [exact, for general a, b; rank drops where a = 0 and b = 0; displayed formula undefined where a = 0]
    solve a*x + y = 1, x + a*y = 0 for x, y
        x = a/(a^2 - 1), y = -1/(a^2 - 1)
        [exact, for general a; rank drops where a^2 - 1 = 0; displayed formula undefined where a^2 - 1 = 0]

Rank drops where **all** nonzero minors of the generic-rank size vanish, on the input domain. A chosen basis
can have extra poles without a rank drop. Substitute particular values and run again when an exceptional case
is needed. Generic conditions survive storage and arithmetic; no random test is labelled a proof.

Scope of this version:
- Fractional exponents in a finite quotient are refused. Surds or i in a denominator continue through the old
  inverse only when that denominator has no ordinary letters. Surds mixed with letters in a denominator are
  refused; a letter denominator requires rational coefficients in the numerator too.
- Explicit surd-series requests such as `1/(sqrt(2)+x) to x^3` still work. Series coefficients retain the existing
  C representation; expansions requiring inversion of a coefficient such as `a+b` remain unsupported.
- Special parameter cases are not solved automatically. Listing all maximal minors can be expensive for large
  rectangular matrices. The shared UP representation retains its existing degree limit of 10000 per letter.
- See [the slice 8 verification record](docs/SLICE8_VERIFICATION.md) for its acceptance key, check counts,
  measured times and the single old fixture line changed in that slice.


## Rational integrals (slice 9)

A rational integral returns a finite area: a rational function plus logarithmic and circular areas.
Every answer is differentiated back exactly before it is accepted. Logs use absolute values; the arbitrary
constant is omitted. Ask explicitly for a series to retain expansion at zero.

    integral(1/(x^2 + 1), x)                 atan(x)  [exact]
    integral(x/(x^2 + 1), x)                 (1/2)*log|x^2 + 1|  [exact]
    integral(1/(x - 1)^2, x)                 -1/(x - 1)  [exact]
    integral(x^2/(x + a), x)                 -ax + x^2/2 + (a^2)*log|a + x|  [exact]
    integral(1/(1 + x), x) to x^4            x - x^2/2 + x^3/3 - x^4/4 + O(x^5)
    apart(1/(x*(x + 1)^2), x)                a sum of three rational fractions, added back exactly

A stored area supports its derivative and substitution of a real number. Partial fractions support the same
operations as their rational sum. General arithmetic on an area is not yet supported.

    let F = integral(1/(x^2 + 1), x)
    d/dx F                                  1/(x^2 + 1)  [exact]
    F(1)                                    pi/4  [exact]
    integral(1/(1 + x), x, 0, 1)             log(2)  [exact]
    integral(1/(1 + x^2), x, 0, 1)           pi/4  [exact]
    integral(1/(1 + x), x, 0, 1) to 60 places
        0.693147180559945309417232121458176568075500134360255254120680  [bounded: 60 places guaranteed]

Definite integrals check the entire closed interval for poles, including the endpoints and poles of even
multiplicity. Reversed and supported algebraic limits work. `use prelude` is not required: named conic
constants use the prelude equations and the existing certified term-rule evaluator with smaller arguments.

Scope:
- After Hermite reduction, the remaining denominator must split into supported linear and quadratic factors.
  An unsupported higher-degree factor is refused with `later: Rothstein-Trager`; bounded factor searches do
  not claim irreducibility. The denominator degree limit is 64.
- Parameter-dependent log/atan parts are restricted to a single linear factor, including its repetitions.
  A quadratic with parameters is refused with a reason. Constant parameter denominators and general linear
  coefficients such as `(a+b)*x+1` use rational coefficients in those parameters.
- Exact definite values may retain named `atan(k)` values beyond the familiar rational multiples of pi.
  General symbolic identities between logarithms, conic-area arithmetic, integration of conic areas, and
  cross derivatives in other letters remain outside this version.
- For algebraic bounds extremely close to a pole, inability to separate them at the working precision is
  reported rather than guessed. As elsewhere, cancelled holes in a rational expression are not retained.

See [the slice 9 verification report](docs/SLICE9_VERIFICATION.md) for all supplied answers, check counts,
changed old outputs, measured test time and remaining limits.

## Notation (slice 10)

What is printed reads back as the same value. Letters are single, as in Newton: `2ax` is 2·a·x unless `ax` is a
defined name. Terms follow the powers of the flowing letter: `x^2 + 2ax + a^2`, `(x + a)/(x - a)`.

    (x^2 + 2ax + a^2)/(x^2 - a^2)     (x + a)/(x - a)  [exact]
    integral(1/(x^4 - 1), x)          -log|x + 1|/4 + log|x - 1|/4 - atan(x)/2  [exact]
    apart(1/(x(x + 1)^2), x)          -1/(x + 1) - 1/(x + 1)^2 + 1/x  [exact]
    integral(1/(x^2 + x + 1), x, 0, 1) - sqrt(3)pi/9     0  [exact]

`tests/roundtrip.sh` checks this for every exact result in the test files.

## Areas under square roots (slice 11)

A letter is put for the root, as Newton did; the root of a quadratic leaves the area of the hyperbola (a
logarithm) or of the circle (an arc):

    integral(sqrt(x^2 + 1), x)          x*sqrt(x^2 + 1)/2 + log|x + sqrt(x^2 + 1)|/2  [exact]
    integral(1/sqrt(1 - x^2), x)        asin(x)  [exact]
    integral(1/(x sqrt(x + 1)), x)      -log|sqrt(x + 1) + 1| + log|sqrt(x + 1) - 1|  [exact]
    integral(sqrt(1 - x^2), x, 0, 1)    pi/4  [exact]
    integral(1/sqrt(x^2 + 1), x, 0, 1)  log(sqrt(2) + 1)  [exact]

A root in a denominator with other factors goes through Euler's substitution (slice 12): a rational point of the
conic s^2 = q gives the new letter t = (sqrt(q) - c)/(x - k), and the area becomes a rational one:

    integral(1/(x sqrt(x^2 + 1)), x)       log|(sqrt(x^2 + 1) - 1)/x|  [exact]
    integral(1/(x^2 sqrt(x^2 + 1)), x)     -sqrt(x^2 + 1)/x  [exact]

Not yet: several different roots, a root of degree 3 or more (elliptic; the series still works), conics with no
rational point found, and at most 32 letters and named constants in one session.

