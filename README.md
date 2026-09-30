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
fractions or negative (a^(1/2), x^(-2/5)). Division is by a single term. `to x^8` names the series letter when
it is not clear.

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

Not yet:
- irrational numbers in coefficients (the starts ±sqrt(2ax) are listed but not followed);
- division by a sum of letters;
- arithmetic on fractional-power results;
- values at approximate arguments.
