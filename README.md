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

Not yet:
- more than one letter in a series;
- Newton's parallelogram for multiple roots;
- fractional powers of x;
- irrational coefficients;
- values at approximate arguments.
