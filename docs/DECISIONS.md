# newtonmath: decisions log
Each decision gets: the question, the generic answer, what Newton did (sources), a recommendation, then the
accepted decision. Mert currently delegates new design decisions to the previous agent: stop and give him a
prompt to relay, then wait for the relayed decision before building. Already approved work may be completed.

## Order of decisions (proposed)
0. Purpose and users (Mert's; the brief leaves it blank).
1. The foundation: the smallest set of primitive objects and operations.
2. Exact vs approximate: how a result states its precision.
3. Checking inside the language: the second route, the verdict.
4. Meaning, then notation: the syntax.
5. The host for the first implementation (bootstrapping towards "no other language").
6. The first slice to build and its tests.
7. Adoption: what makes others use it (Leibniz's notation outlived Newton's).

## D1 research: the foundation
- NATP00180 (1676): he found the series by interpolating Wallis's sequences (pattern), then "neglexi penitus
  interpolationem serierum et has operationes tanquam fundamenta magis genuina solummodo adhibui".
  He dropped the pattern and kept division and root extraction as the foundation.
- NATP00180: "vix speraverim dari posse aut simpliciora aut magis generalia fundamenta ... quam sunt divisiones,
  et extractiones radicum".
- NATP00180: "res non eget demonstratione prout ego operor. Habito meo fundamento nemo potuit tangentes aliter
  ducere, nisi volens de recta viâ deviaret". With the right foundation, the method needs no separate proof.
- NATP00182 / NATP00296: operations on letters are done as in decimal numbers ("in speciebus ac ... in decimalibus").
- NATP00204, rules I-III: primitive, composition, reduction to simple terms.
- NATP00296: the two problems, direct (fluents -> fluxions) and inverse ("et vice versa").
- M12 (De analysi): infinite equations are as legitimate and exact as finite ones.

## D0 decided (Mert)
- Users: mathematicians who want to do mathematical computation.
- Scope: arithmetic, algebra, analysis. Arithmetic geometry yes. Simple geometry if it comes easily.
  Logic only where it helps.
- Program it the way Newton would program; consult Newton at every decision; read in depth.
## D1 accepted with D0's scope
- One primitive object, number, with its extensions: rational, letter, series, flowing quantity.
- Primitive operations: the four operations, root extraction and resolution.
- The direct/inverse pair; everything else from rules I-III.

## D2 research: exact vs approximate
- NATP00296 (Methodus, Problem 1): substitute x+mo; terms with o "respectu cæterorum nihil valebunt", rejected.
  M47 (1704): later he rejected this: errors, however small, are not to be neglected; prime and ultimate ratios.
  Lesson: working modulo o^2 is exact; neglect is never silent.
- NATP00204 (De analysi): the start has a stated bound ("2 numerus qui minùs quàm decimâ sui parte differt a
  radice"). The digits gained are fixed by a rule ("donec tot eliciantur figuræ quot locis ... distant"). The next
  pass gives "bis tot figuras ... unâ dempta". Approximate values are marked "ferè".
- NATP00204: when unsure whether the linear step suffices, keep the quadratic term; explicit criterion (the square
  of the penultimate coefficient at least ten times the product of the last and antepenultimate).
- NATP00204: "operationem sic produco quosque placuerit": precision on demand.
- NATP00204: the residual equation stays, "unde examinatio operis hic æque poterit institui ac in reliqua
  Arithmetica".
- NATP00180: logarithms. A few values (0.8, 0.9, 1.1, 1.2) from fast series; log 2, log 10 and the primes
  7, 13, 17, 37 follow by exact relations. "Pudet dicere ad quot figurarum loca ... perduxi": beyond need is waste.
- NATP00182: roots can be extracted "modis infinitis ... de quorum simplicioribus ... semper consulendum est".
- M19: series modelled on decimal fractions.
Synthesis: decimals (base 10), power series (base x) and p-adics (base p) are one kind of object, a number known
to K places. Exact objects are exact; an approximate one carries its guaranteed places, its residual, and the
rule that fixed them.

## D2 decided (Mert): accepted as proposed
- Exact by default.
- No bare floats: an approximate value carries its guaranteed places, its rule and its residual.
- Decimal, series and p-adic are one "number known to K places".
- Precision on demand.
- Approximation steps test their own safety criterion.

## D3 research: checking inside the language
- NATP00180: the binomial series, found by pattern, checked by the inverse operation: "ut probarem has operationes
  multiplicavi ... in se, et factum est 1−xx ... Quod ut certa fuit harum conclusionum demonstratio, sic me manu
  duxit ad tentandum e converso". The check was a demonstration, and it led to the second route (extraction).
- NATP00296, Problem 2: the inverse is checked by the direct: "si ... potest regredi ... rectè operatus est; sin
  secùs, vitiosè". His counterexample shows a wrong inverse caught this way.
- NATP00204: the residual equation, "examinatio operis hic æque poterit institui ac in reliqua Arithmetica".
  Checking costs as little as in arithmetic.
- NATP00204, end: a general proof for the method (the last term becomes less than any given quantity, by
  Elements 10.1; a geometric bound when x < 1/2), then the case check: "Idem patebit substituendo quotientem pro y
  in æquationem propositam".
- NATP00204: "Demonstratio ejus ex ipso modo operandi patet": the proof is visible in how the method works.
- M35: analysis in private, presented in the most certain (synthetic) form in public.

## D4 research: notation
- NATP00182: "sicut Analystæ pro aa, aaa ... scribere solent a^2, a^3, sic ego pro sqrt(a) ... scribo a^(1/2) ...
  & pro 1/a ... a^(−1)". He extended the notation readers already had, by analogy.
- NATP00182: "pro terminis inter operandum inventis in quoto, usurpo A, B, C, D": names for intermediate results;
  each term is defined from the previous one (B = (m/n)·A·Q).
- NATP00296: roles by letter: a, b, c known; v, x, y, z flowing; l, m, n, r their velocities; o the moment.
- NATP00204: the "maximè expeditum" form for substitution, (((y−4)×y:+5)×y:−12)×y:+17, with ":" closing a group.
  Notation chosen for the computation.
- NATP00101: "Ne hujusmodi operationes obscuræ nimis evadant, Lemmata 6 sequentia ... præmittam". State what the
  work relies on before it starts.
- NATP00204 / NATP00182: computations laid out as tables ("Diagrammata ubi dextra columna prodit substituendo in
  media columnâ").
- M40 and texts/INDEX: his dots (1691) lost to Leibniz's d/dx.
  General knowledge, not from these texts: d/dx names the variable, and rules like the chain rule read as algebra
  on the symbol.

## D3 decided (Mert): accepted as proposed
## D4 decided (Mert): accepted, with an amendment: everything must be typable on a plain keyboard
Research for the amendment:
- NATP00296 / NATP00182: he replaced the radical sign with vinculum (two-dimensional) by linear exponents:
  "usurpo ... x^(1/2) ... x^(1/3) ... pro sqrt()x ... sqrt()c:x"; a^2·(a^3+bbx)^(−1/3) instead of a fraction under a
  cube root.
- NATP00204: ":" closes a group inline ((((y−4)×y:+5)×y:−12)×y:+17) instead of an overbar.
Rule:
- The canonical form is plain ASCII: x^(1/2), sqrt(x), integral(f, x), d/dx.
- Unicode forms (√, ∫, ², ≤, →) are optional aliases that the reader accepts.
- Output is ASCII by default.

## D5 research: the host of the first implementation, and how the user reaches the language
- H2 (constitution): he makes the tool when none is adequate. Instruments letter (NATP00003): he devised his own
  metal polishing ("a tender way of polishing proper for mettall") and made the telescope himself.
- Collaboration line: hand execution to better hands and instruments, keep the judgement; learn the operation
  before directing it.
- NATP00180: when Mercator's Logarithmotechnia appeared, "cœpi ea minùs curare". He did not redo work that others
  had already done well.
- NATP00204: "Quicquid laboris hic est in substituendo ... reperietur": he found where the labour lies and built
  the fastest form for exactly that part (the nested substitution).
Available here: gcc 13, clang, make, rustc, go.

## D5 decided (Mert): B
- The interpreter is written in C (the lathe).
- The number core is our own, with no dependencies (the mirror).
- Linking GMP later would be a separate decision, taken only on a measured bottleneck.
- Self-hosting: the library first, the interpreter later if ever.
- Code lives in the scratchpad until Mert decides on a repository.

## D6 research: the first slice and its tests
- NATP00128 <2r>, the first page of his notebook: "Of the extraction of Pure Square Cubick ... rootes". Digit-group
  pointing, the "first side" A, the "second side" E, the divisor. Worked example: square root of 2916 = 54.
  His own order starts from exact integer arithmetic and root extraction.
- NATP00204: the resolution is shown first in numbers (y^3−2y−5=0 → 2.09455147), then "His in numeris sic
  ostensis: Sit æquatio literalis ...". Numbers first, letters second.
- NATP00296: the Methodus order: reduction by division → root extraction → affected equations (numbers, then
  species) → the problems.
- Constitution, new-problem move: do the simplest instance exactly, then generalise.
- Known answers in the texts, usable as tests:
  - sqrt(2916) = 54 (NATP00128);
  - y^3−2y−5=0 near 2: 2.09455147 (NATP00204; true 2.0945514815...; his last printed digit is off by one);
  - 1/(1+x^2) = 1−x^2+x^4−...; sqrt(aa+xx) = a + xx/2a − x^4/8a^3 + ... (NATP00296);
  - (1−xx)^(1/2) = 1 − x^2/2 − x^4/8 − x^6/16 − 5x^8/128; (1−xx)^(1/3) = 1 − x^2/3 − x^4/9 − 5x^6/81 (NATP00180);
    each checked by raising back to 1−xx;
  - log 2 = 0.6931471805597, log 10 = 2.3025850929933, from the logs of 0.8, 0.9, 1.1, 1.2 (NATP00180);
  - the literal root of y^3 + aay − 2a^3 + axy − x^3 = 0: y = a − x/4 + xx/64a + ... (NATP00204).

## D5b decided (Mert): (a)
C first. Later newtonmath gets its own compiler that emits machine code (self-hosting).
## D6 decided (Mert): slice 1 = numbers
- Exact integers and rationals.
- Root extraction to guaranteed places.
- Numerical resolution of equations.
- Parser, printer, REPL and file mode.
- Tests from Newton's worked examples plus inverse checks, ms each, under 60 s in all.

## D7 research: the internal form of a big number
- NATP00128 <2r>: "Let the number whose roote is to bee extracted bee pointed ... comprizeing soe many numbers
  under each point as the number hath dimensions": decimal digits grouped into blocks (2 for square roots, 3 for
  cubes, 5 for fifth roots). Example 5708·63524·10802 (fifth root) and 29·16 → 54.
- M19 / NATP00296: series modelled on decimal fractions.
- D2: "a number known to K places" must truncate at exact decimal places.
Options:
- A: blocks of 9 decimal digits (base 10^9): Newton's pointing; exact decimal places; printing is trivial.
- B: binary words (base 2^32 or 2^64), as GMP does: multiplication about 2-3x faster; printing needs a conversion;
  decimal places need extra work.

## D7 revised after Mert rejected "redo at double precision" (not a good architecture)
Research: Newton fixes the precision BEFORE and DURING the work; he never repeats it.
- NATP00296: "circa finem operis eos omnes terminos negligo quorum dimensiones transcenderent dimensiones ultimi
  termini ad quem cupio quotientem solummodò produci". The target is chosen first. Every intermediate step is
  carried only as far as the target needs.
- NATP00204: "dividendo ... donec tot eliciantur figuræ quot locis primæ figuræ hujus & principalis quotientis
  exclusivè distant". The number of certain figures follows from the positions of the leading figures, which are
  known cheaply beforehand.
- NATP00204: the next pass gives "bis tot figuras ... unâ dempta". The loss is counted in advance.
- NATP00204: "operam sub fine ... levabis, si primos ejus terminos gradatim omiseris". Each step is done at the
  precision that step can deliver: little at the start, more as the digits grow.
- NATP00204: the quadratic-term criterion, decided from coefficients already at hand before the step is taken.
- NATP00204, end: a tail bound (x < 1/2 ⇒ each term exceeds half the rest), so the number of terms is known
  before summing.
- NATP00180: logs from the fastest series (arguments 0.1, 0.2), the rest by exact relations.
Architecture:
- A plan pass: from the requested places, go backwards through the expression. Each operation has a known
  digit-loss rule and a magnitude, obtained from cheap leading digits, so the working places are assigned
  in advance.
- One compute pass, where each Newton step runs at its own growing precision (the "gradatim" rule).
- One final check of the result's radius and residual. A failure is a verdict (a bug or an unforeseen
  cancellation), never a silent retry.
- Where magnitudes cannot be known in advance, e.g. a difference that may be exactly 0: keep the value exact
  (algebraic) and approximate only at output.

## D7 decided (Mert)
- Base 10^9 blocks.
- Exact by default; a root is kept as its equation plus a verified bracket.
- Approximate values use plain ball arithmetic with fixed guard digits and local radius rules.
- Roots and series continue from their stored state, never recompute.
- An expression that cannot reach the requested places says how many are guaranteed; no silent retries.
Surveyed: PARI (fixed precision), Arb (Ziv retries), Mathematica (significance arithmetic, $MaxExtraPrecision),
SymPy evalf (adaptive), Calcium (exact fields with enclosures).

## Slice 1 built (scratchpad/newtonmath)
C, no dependencies, 1300 lines. `make test`: 42 known-answer checks plus 9065 integer-core checks, about 1 s.
Small choices made while coding (reported to Mert):
- default 20 places;
- 30 guard digits for balls;
- side-by-side multiplication (2y);
- an unbound letter is an unknown (a polynomial);
- a sign-change bracket widens up to ±500 units in the last working place;
- the real cube root of a negative number.

# Slice 2: letters and series. Research
- NATP00296 (opening): series with x^(1/2), x^(3/2), x^(−1) from the start, "ob analogiam rei". The ascending
  order can be inverted: "Potest etiam ordo terminorum inverti ... xx+aa, et radix est x + aa/2x − a^4/8x^3".
  The series object has rational exponents and can be taken in x or in 1/x.
- NATP00296: "circa finem operis eos omnes terminos negligo quorum dimensiones transcenderent ... ultimi termini ad
  quem cupio quotientem ... produci". Truncation is named by the last degree kept.
- NATP00296 and NATP00180 (Epistola posterior): the parallelogram. For y^6 − 5xy^5 + (x^3/a)y^4 − 7a^2x^2y^2 +
  6a^3x^3 + b^2x^4 = 0, the ruler rests on x^3, x^2y^2, y^6. From y^6 − 7a^2x^2y^2 + 6a^3x^3 he gets four starts
  ±sqrt(ax), ±sqrt(2ax), after the reduction v^6 − 7v^2 + 6 = 0 with y = v·sqrt(ax). "quorum quemlibet pro
  primo termino Quotientis accipere [liceat]": the user chooses the branch. "y denotat radicem extrahendam et x
  alteram indefinitam quantitatem ex cujus potestatibus series constituenda est": roles are named.
- NATP00204: the literal resolution y^3 + aay + axy − 2a^3 − x^3 = 0 → y = a − x/4 + xx/64a + ...; checked by
  substitution: "Idem patebit substituendo quotientem pro y in æquationem propositam".
- NATP00180: "Denominatores fractionum ... reducantur ad quam paucissimas et minimè compositas". Coefficients
  with simple denominators.
- NATP00182: the binomial rule with named terms A, B, C, D, each from the previous (B = (m/n)·A·Q). The term
  ratio is known, so the tail can be bounded.
- NATP00204, end: the tail bound (x < 1/2 ⇒ each term exceeds the rest), and "numeros coefficientes ... plerumque
  decrescent": the bound needs the coefficients' behaviour.
- NATP00180: logs of 0.8, 0.9, 1.1, 1.2 from the area series, then log 2 and log 10 by exact relations.
- De analysi rules I-II: area (integral) term by term.
- M14: reversion of a series (the sine from the arcsine).

## Slice 2: the modern world consulted (Mert: Newton must also consult the modern world)
- Series objects:
  - PARI t_SER: integer exponents, precision O(x^n).
  - Sage PuiseuxSeriesRing: finite precision, rational exponents.
  - Sage LazyPowerSeriesRing: terms computed on demand. This is Newton's "quousque placuerit" and our D7
    continuation.
  - Maple algcurves[puiseux] and Mathematica also give Puiseux series.
- Branches: Maple's puiseux returns all branches; Newton lets the user choose one.
- Coefficients with letters: FLINT fmpz_mpoly_q keeps full multivariate rational functions canonical by gcd.
  An open FLINT discussion (issue 1774) proposes factored denominators to avoid a gcd at every operation. This is
  close to Newton's "denominators as few and least compound as possible".
- Tail bounds:
  - SymPy hypsum and Arb use the known term ratio (hypergeometric series).
  - Mezzarobba (NumGfun; Ann. Henri Lebesgue 2019) gives majorant-series truncation bounds for D-finite series,
    the general case.

## Slice 2: what Newton would do with the modern capabilities
- NATP00296, Problem 2, Example 1: the fluxional equation n/m = 1 − 3x + y + xx + xy (y' = 1 − 3x + y + x^2 + xy,
  y(0) = 0), solved in a table term by term: y = x − xx + x^3/3 − x^4/6 + x^5/30 + ... "Et eâdem operatione
  sæpiùs repetitâ Quotientem ad arbitrium producere possis": on demand, as far as wanted.
- The modern D-finite theory (linear ODEs with polynomial coefficients, recurrences, Mezzarobba's bounds) is his
  Problem 2 made general and rigorous. His binomial rule (A, B, C, D) is a P-recursive recurrence.
- Synthesis, one consonant frame:
  - a number is known by its equation (slice 1: a root is its polynomial plus a bracket);
  - a function is known by its fluxional equation plus its start.
  From that one object come terms on demand (lazy), tail bounds, values at points, the check by substitution, and
  closure under +, ×, d/dx, integral. Algebraic functions (Puiseux branches) are D-finite too.
- Nonlinear equations stay with Newton's resolution by series, with an honest verdict (no general tail bound).

## Slice 2 redesigned by Newton's behaviour (Mert: no over-engineering, no guardrail layers; do it right from the
## start; be bold; copy his behaviour, don't get lost in modern ideas)
His behaviour, from the texts:
- He reduced all of analysis to one object and a few operations: infinite series, handled "in speciebus ac in
  decimalibus numeris". A decimal is a series in powers of 1/10, a function a series in powers of x.
  Same operations, same work.
- He named the whole method in one sentence: given an equation with fluents, find the fluxions, and conversely.
- He then used one engine for everything hard: resolution. Start from a known value, substitute, keep the
  residual equation, and gain places as far as wanted. The check is the residual itself, not a separate layer.
- He did the simplest case exactly and immediately stated the general rule (the binomial at once for every m/n).
Design:
- One engine, resolution, over two domains: decimal places and series degrees.
- Verbs: + − × ÷, ^ (fractional powers are root extraction), `root of <equation>` for algebraic and fluxional
  equations alike, d/dx and integral.
- Named functions (log, exp, sin, ...) are not written in C. They are defined in newtonmath by their equations.
- Output: the value, how far it is carried, and what it satisfies.
Deferred until a problem needs them: letters with factored denominators, Puiseux branches (the parallelogram),
general tail-bound theory.

## Slice 2 built
- One resolution engine, used for numbers (Newton's iteration on places) and for series (term by term with the
  moment o, o·o rejected).
- `root of` serves algebraic and fluxional equations. d/dx and integral are exact on polynomials and term by term
  on series.
- Fractional powers use the rule that generalises the binomial, checked by raising back.
- Named functions live in lib/prelude.nm, written in newtonmath.
- Found by the engine: Newton's Methodus Problem 2 table gives +x^6/45; substitution gives −x^6/45. By hand,
  6c6 = c5 + c4 = 1/30 − 1/6, so c6 = −1/45.
Small choices made while coding:
- default order x^8;
- `to x^N` keeps degrees up to N;
- a series letter is the one free letter, x if none;
- a defined series is stored as its definition and re-derived at the order asked.

# Slice 3: the value of a series at a number. Research (Newton's behaviour)
- NATP00180, the logarithms:
  - He computes only at small arguments (±0.1, ±0.2), where each term is about 100 times smaller than the last.
    He uses the half-sum and half-difference, x + x^3/3 + x^5/5 + ... and x^2/2 + x^4/4 + ..., so that only
    every other term is needed.
  - Everything else follows from exact relations: log 2 from 1.2·1.2/(0.8·0.9) = 2, log 10 from 2·2·2/0.8 = 10,
    and the primes 7, 13, 17, 37 from 0.98, 0.99, 1.01, 1.02.
  - He prints log 2 = 0.6931471805597 (true 0.69314718055994...) and log 10 = 2.3025850929933
    (true 2.30258509299404...).
- NATP00182 / NATP00180: each term from the one before (A, B, C, D). He knows the term rule, so he knows the tail.
- NATP00204, end: the tail is bounded once the terms shrink geometrically (x < 1/2 ⇒ each term exceeds all the
  rest together).
- NATP00180: "Pudet dicere ad quot figurarum loca ...": carry the work only as far as needed.
Design:
- The value of f at a rational number comes from f's term rule. The rule is derived mechanically from f's
  equation when the equation is linear in y, y', y'', ... with polynomial coefficients: Newton's Problem 2 table
  written as a rule.
- The terms are summed exactly, and the tail is bounded from the rule: with |coefficient ratio| ≤ β_t for k ≥ K,
  the tail is at most s·γ·W/(1 − γ), where γ = Σ β_t r^t < 1. The sum stops as soon as the bound is below the
  places asked. There is one pass and no retries.
- The choice of argument stays with the user, as with Newton. Near the edge of convergence (γ ≥ 1) the answer
  is a plain statement that the terms shrink too slowly there.
- The prelude is rewritten so that every function is defined by such an equation:
  - exp: y' = y
  - sin, cos: y'' = −y
  - log1p: (1+x) y' = 1
  - atan: (1+x^2) y' = 1
  - asin: (1−x^2) y'' = x y'

## Slice 3 built
- The term rule is read off the equation. y and its derivatives are set to 0, 1, 2 to get the polynomial
  coefficients, with an exact linearity check, and turned into D(n) c_n = sum N_t(n) c_{n-t} − h_{n−s}.
  The rule must reproduce the resolution's own terms before it is used.
- Tail bound: |N_t(n)/D(n)| ≤ U_t(N) for n ≥ N, where U_t is non-increasing, built from absolute values of
  coefficients. The tail is ≤ T·γ·W/(1−γ) with γ = Σ U_t(N) r^t < 1. If lim γ ≥ 1, the answer is that the
  terms shrink too slowly.
- The terms are summed as balls (exact coefficients, powers of the argument carried as places), in one pass.
- A let-bound approximate number is kept as its recipe and carried again to the places later asked.
- The prelude is now all linear equations: exp, sin, cos, log1p, atan, asin.
Timing:
- log 2 by Newton's relations: 100 places in 10 ms, 1000 places in 1.8 s.
- e: 1000 places in 0.14 s.
Tests: log 2, log 10, e, sin, cos, atan, pi (twice) and exp(−3/2), each to 30 and 200 places against bc -l.

## Slice 4 built
- Quantities in letters: sums of k·a^e1·b^e2 with rational k and rational exponents. Division by a single term
  only. Polynomials in the exact pass are now such quantities (several letters, Laurent: a^2/x is exact).
- The series letter follows Newton's convention when several letters appear: initial letters are given, final
  letters (t, u, v, w, x, z) flow; otherwise `to x^N` names it.
- The parallelogram (src/newton.c):
  - the lower edges of the marks (j, i) give the starts y = c x^g;
  - the ruler's equation is reduced to numbers by Newton's scaling (y = v·M) and solved for rational v;
  - other roots are listed with their equation, marked as later;
  - with no start, the starts are listed;
  - a simple start continues term by term, with z(0) = 0 and each coefficient −F_d / F_z(0,0);
  - a multiple start applies the parallelogram again;
  - the finished root is substituted back into F in u = x^(1/D).
Found in Newton's text: for x^2y^5 − 3c^4xy^2 − c^5x^2 + c^7 = 0 he writes the start as the fifth root of
c^7/xx. The terms on the ruler, x^2y^5 + c^7, give y^5 = −c^7/x^2, so the real start is −c^(7/5)x^(−2/5).
Matches Newton: the literal resolution (all five printed terms); the start 3x; the four starts of the y^6
example (±sqrt(ax) followed, ±sqrt(2ax) listed as irrational).

# Slice 5: toward the language's own compiler, by the mathematical road
Decision (Mert): (a). The language stays for mathematics. It gains Newton's way of writing repetition, so that
its methods can move from C into the language one by one. A small kernel (numbers, reader, printer) stays in C
until a compiler for it comes last. Newton: the foundation is the smallest set of operations (M06); repetition
is a rule ("eâdem operatione sæpiùs repetitâ"; each term from the one before, A, B, C, D); tables are built
forward (M26); only the bottleneck part of the tool is made by hand.
Built:
- `let f(x, y) = ...` rules, with lexical frames;
- `let A[0] = ..., A[k] = ...` sequences, with tables built forward and invalidated by any new definition;
- `v1 if c1, ..., v otherwise` cases;
- `sum(... for k = a to b)`.
Proof of the road: Newton's binomial rule, his resolution of y^3 − 2y − 5 and his Problem 2 table are written in
the language, and each agrees with the engine built into C.
Small choices made while coding:
- rules may call themselves up to 4000 deep; sequences have no such limit;
- the stack limit is raised to 256 MB at start;
- a sequence's terms must be exact.

# Slice 6: surds and i
Newton (Methodus, Problem 1, Ex. 3): "Siquando in æquatione propositâ insint fractiones complexæ aut surdæ
quantitates, pro singulis pono totidem literas ... operor ut ante. Dein supprimo et extermino literas
ascriptitias." In the parallelogram he writes y = v sqrt(ax) with v^6 − 7v^2 + 6 = 0. In the notebook, sqrt(−b)
is "impossible"; here it is i.
Built:
- a surd is a letter with its equation and a chosen root:
  - a certified bracket, for roots of rational numbers and named roots;
  - a radical of a positive quantity, for surds over surds;
  - i.
- Products are reduced by the equations.
- Inverses: a single term directly; otherwise Euclid between the quantity and the latest surd's equation, with
  coefficients in the earlier surds. If the equation factors, the factor the root satisfies is kept.
- Radicals of integers take out small n-th powers: sqrt(8) = 2sqrt(2).
- Values: balls from the brackets. A surd over a surd takes the roots of both ends of its radicand's ball.
  Complex quantities print as a + bi.
- The value engine accepts approximate arguments (an upper bound of |a| for the tail).
- The parallelogram accepts a start with surds: it is substituted into the ruler's equation.
Small choices made while coding:
- a bare `i` is the square root of −1 (a sum index or rule parameter named i still works);
- a lone surd prints as itself, `sqrt(2) [exact]`, and gives its places when asked;
- sqrt(2)·sqrt(3) stays a product of two surds (not sqrt(6)).

# Slice 7: four modules, then one integration pass
Decision (Mert): build systems, linear algebra, complex arguments and number theory each apart, then join them to
the engine once. Newton was consulted instead of Mert at each choice, taking a practical method where it fits his
way better than his own historical one.
Newton: from two equations one unknown is exterminated (Arithmetica Universalis, "De duarum pluriumve
aequationum in unam transformatione"); every root found is put back into the equation; a moved series
(Methodus: y at a + x) has its starting values as given quantities, written as letters.
Built:
- `solve eqs for x, y`: extermination by resultants, roots by rational roots, the quadratic formula and Sturm
  brackets (r_1, r_2 ... with their equations and certified places), every solution put back;
- `eliminate y from A = B, C = D`: the resultant, which can be named with let;
- matrices: literals, + − ×, scalar × and /, whole powers (negative through the inverse), det, inverse,
  transpose, rank, nullspace, charpoly, eigenvalues, linsolve (with the null space when not unique);
- whole numbers: gcd, lcm, mod, powmod, invmod, isprime, factor, divisors, sigma, phi, nextprime;
- values at complex points and f'(a), f''(a) ... from the term rule (one routine, rule_value, for all);
- f(a + x) as a series: an equation-defined f is moved to a, with f(a), f'(a) ... as named letters; an
  expression-defined recipe has the series (or number) put in.
Small choices made while coding:
- built-in names are used only when the name is not defined by the user;
- the limit on different letters and surds in one session went from 12 to 32 (a session of matrices with
  letters, roots and surds ran out at 12; the tests run about 25% slower);
- a solve prints its own legend for r_k, and none for radicals like sqrt(5), which explain themselves;
- `let s2 = sqrt(1 + x)` then `s2(1/2)` now gives sqrt(6)/2 exactly (before: an error asking for an equation).

## Evaluation and finite-sum order
Question: should compiler argument order or earlier statements decide how an expression is evaluated and
how equal-degree terms are printed?
The generic answer: C leaves function-argument evaluation order unspecified. A language must choose its own
evaluation order; a printer can order terms independently of the internal representation.
Newton:
- NATP00182: the root computations are laid out as tables, substituting the values from the left column into
  the middle column to obtain the right column. The computation has explicit steps.
- NATP00296: after extracting the root of aa + xx, he explicitly reverses the term order to xx + aa and
  expands about the other leading term. His ordering serves the calculation; it does not prescribe a universal
  alphabetical printer or a programming-language evaluation rule.
Recommendation: evaluate operands from left to right; print finite sums by descending total degree, then
descending exponents in ascending alphabetical letter-name order. These are practical implementation choices.
Decision (Mert): approved both rules.
Built:
- binary expressions evaluate their left operand before their right operand;
- an undefined name applied to an argument evaluates its named factor before the argument, in both the exact
  and series passes;
- finite-sum printing uses letter names to break degree ties, rather than registration indices;
- regression checks cover both registration orders, full names, fractional/negative exponents, and the first
  error in a binary expression.
Small choices made while coding:
- build the printer's alphabetical index order once per sum; compare exact rational degrees and exponents;
- keep the arithmetic's indices and surd equations intact; the printer sorts a copy of the terms;
- review only the two changed fixture lines: the binomial square and the resultant obtained by substituting
  y = x + b into y^2 = x^3 + a. Both were checked independently by substitution and with bc.

# Slice 8: rational functions

Question: how should division by a sum of ordinary letters behave, and how can matrices use it?
Decision (Mert, supplied design D1-D6): return a finite exact quotient first; expand only when a series is
requested. Reduce by recursive multivariate polynomial gcd over Q, using primitive pseudo-remainder sequences,
not modular or heuristic gcd. Share elimination's univariate-in-one-letter representation. Matrices use the
resulting field; rank and nullspace must state their generic meaning and exceptional values. Check quotients by
cross-multiplication and every gcd by exact division of both inputs. The supplied design supersedes the earlier
requirement to return new design questions to the previous agent; routine implementation choices are delegated.

Newton and evidence:
- NATP00204, De Analysi, rule III and the following division example (local text, lines 27-33): aa/(b+x) is a
  given quantity, expanded by division to compute its area. This motivates keeping the finite form available.
- NATP00296, Problem 2: checking the inverse by the direct operation, with a failed check called "vitiose"
  (the passage already recorded under D3 above).
- Arithmetica Universalis: fraction reduction by continual division / Euclid on polynomials, **from general
  knowledge**; this book is not in sources/ and this historical attribution has not been checked against it.

Isolated module built:
- `R { C num, den }`, arithmetic, whole powers, equality, persistence and printing in src/ratfun.c.
- Recursive content gcd and primitive pseudo-remainders use UP allocated, split and subtracted by elim.c.
- Exact polynomial division multiplies back. Each gcd is checked against both inputs; quotient normalization
  checks the original numerator and denominator by cross-multiplication. Unit-denominator sums/products use
  the existing C operations directly.
- Unsupported fractional powers in quotients, surds mixed with ordinary letters in denominators, and zero
  denominators fail with a reason. Existing inverses of pure surd/complex denominators still use c_inv.

Small choices:
- Main letter: alphabetically last ordinary letter present. Leading term: the existing fixed descending total
  degree, then descending exponents in alphabetical letter-name order; the printer shares this comparison.
- Make the gcd and denominator monic over Q, stronger than merely positive. This also fixes rational scalar
  units and gives a unique normal form independent of registration order. gcd(0,0)=0; gcd(0,p) is monic(p).
- Clear negative whole powers in numerator and denominator by the same monomial before polynomial gcd.
- A quotient with an ordinary-letter denominator requires rational coefficients; surd numerators in such
  quotients are outside v1. An ordinary-letter numerator over a pure surd denominator can still be rationalized.
- Tests use explicit multiplication in input expressions: `a*x` is a product; `ax` remains a whole name.
- Numerical checks are fixed-seed regression evidence, not proofs: 64 nonsingular rational assignments compare
  11 original expressions evaluated with Q alone against the reduced results (704 comparisons).

Integration built:
- V_RAT participates in exact arithmetic, whole powers, equality, rules, sequences, printing and permanent
  storage. A constant denominator collapses to V_Q/V_POLY. Integer-only gcd keeps its old behavior; a
  polynomial argument selects polynomial gcd. The rest of number theory is unchanged.
- Explicit `to x^N` enters the existing series pass immediately. This preserves the old surd-series example
  `1/(sqrt(2)+x) to x^3` while the corresponding finite quotient is refused by v1's scope.
- A bound rational quantity in one ordinary letter supports exact substitution and series substitution. Its
  numerator and denominator use UP Horner evaluation. Exact derivatives use the quotient rule; integration
  still uses the existing series path. Prime syntax on a bound rational quantity is refused with instructions
  to differentiate first; it is not silently ignored.
- Mat entries and their operations use R. Berkowitz's recurrence and the existing elimination procedure are
  retained; inverse, solutions and null vectors are multiplied back. Whole matrix powers avoid signed overflow
  when taking the magnitude of a negative exponent.
- Parametric affine `solve` statements use mat_solve. Numeric and nonlinear systems retain the elimination
  solver. Inconsistent parametric systems report inconsistency for general values; exceptional solutions are
  not classified.
- Rank-drop conditions are the simultaneous zeros of all nonzero minors of generic-rank size, on the input
  domain. A nonzero constant minor prevents a drop. This is exact and deterministic, though listing all such
  minors can be expensive for large rectangular matrices. A failed chosen pivot is not a rank-drop condition.
- Input-denominator poles and poles of a displayed solution/nullspace basis are printed separately. For [a,b],
  rank drops at a=b=0, while the chosen basis [-b/a,1] has a pole at a=0 even if b is nonzero.
- Generic qualifications persist with bound values and flow through arithmetic, builtins and matrix entries.
  A rounded or bounded output keeps its own verdict alongside these qualifications.
- Equality of rational functions is equality in Q(letters); cancellation does not retain the removed holes
  of an original expression (p/p=1 as required). No classification of all exceptional parameter cases is added.
- One old fixture changes: tests/integration.out line 70, the binding g=1/(1-x), now prints -1/(x-1) [exact]
  instead of its default geometric series through x^8. Multiplying by 1-x gives 1; its explicit expansions
  and its value g(1/3)=3/2 are unchanged. No other pre-existing .out line changes.

# Slice 9: rational integrals and conic areas

Question: how should the finite area of a rational curve be represented, computed and checked?
Decision (Mert, supplied D1-D7): polynomial division, Hermite reduction, then linear/quadratic factors over Q;
return a rational part plus logarithmic and circular areas. Keep explicit series requests unchanged. Add
`apart(f,x)`, exact definite integrals, numerical substitution, and guaranteed decimal values through the
existing equation/term-rule engine. Every primitive is differentiated back and every partial fraction sum is
added back using exact arithmetic. General arithmetic on the new area value is outside this slice.

Newton and evidence:
- The supplied design attributes reduction to conic areas to Methodus, Problem IX and its tables. The local
  NATP00295 text is a 115-line appendix on geometric fluxion axioms, and NATP00296 is a 537-line Methodus
  transcription ending before Problem IX. Those local files do not verify the specific table attribution;
  it is the supplied design's historical motivation, not a newly verified historical claim.
- NATP00296, Problem 2, local line 184 explicitly checks the inverse result by the direct problem and calls a
  disagreement "vitiose". This directly supports the implemented derivative check.
- Hermite reduction and the bounded exact factor search are practical modern implementation choices under
  the approved method; they are not attributed to Newton.

Representation and small implementation choices:
- `Integral` stores its integration letter, rational part R, and terms `(R coefficient, LOG|ATAN, C polynomial)`.
  R coefficients allow `1/((a+b)*x+1)` without introducing reciprocals of sums into C. No integration constant
  is printed. `log|p|` is the real logarithm of the absolute value, on intervals where the primitive is defined.
- Polynomial division comes first. With D=G*S, G=gcd(D,D'), Hermite reduction solves
  `A = S*B' - (S*G'/G)*B + G*C` by exact coefficient comparison using the existing matrix solver.
  The rational primitive is B/G; only C/S is sent to conic factorization. A zero remainder can therefore
  integrate rationally even when D contains an otherwise unsupported cubic.
- Existing rational-root routines are shared with elimination. The sample roots 0, -1, 1 are tried regardless
  of coefficient size. Quadratics use the discriminant. For larger factors, a bounded Kronecker search
  interpolates integer quadratic candidates from divisors of the values at -1, 0, 1 and divides each back.
  It can find products of quadratics with no rational roots. Failure never claims a proof of irreducibility.
- Negative discriminants give atan plus, if needed, a log of the quadratic. Positive nonsquare discriminants
  give conjugate surd logs. Private unreduced derivative fractions allow conjugate cancellation before using
  the ordinary Q-rational normalizer; the general R arithmetic scope is unchanged.
- The denominator degree is limited to 64. Rational-root trial search retains the existing coefficient bound
  below 10^12; quadratic interpolation uses at most 256 divisors per sample and 100000 signed candidates.
  Search exhaustion is reported, with `later: Rothstein-Trager` for an unsupported remaining factor.
- `apart` stores its individual rational summands for display and persistence; ordinary arithmetic,
  differentiation, numeric substitution and explicit series work on their checked sum.
- The new area value supports persistence, derivative in its integration letter, and substitution of a real
  number. Another integral, another derivative letter, or general arithmetic on an area gets an explicit
  refusal. Purely rational primitives collapse to the existing rational/polynomial value types.
- For nonrational inputs, the existing polynomial power rule and series resolver remain in charge. An explicit
  `integral(f,x) to x^N` always enters the old series route. Definite integrals cannot silently discard their
  bounds when given a formal series order.

Definite values and checks:
- Endpoints are exact real numbers, including supported surds. Exact substitution detects endpoint poles;
  Sturm variation counts detect interior poles, including repeated ones. Rational endpoints are tested
  exactly; algebraic endpoints use enclosing rational balls after exact endpoint substitution. An unresolved
  separation near an algebraic endpoint is refused. Reversed intervals keep their orientation.
- Pole checks use the reduced rational function, consistent with slice 8: cancelled holes are not retained.
  Equal endpoints give zero only when that point is regular. Unknown parameter-dependent poles are refused.
- Evaluated logs and atan values are named numbers with callbacks, through the same fluxional equations as
  `lib/prelude.nm` and the existing `series_value_d`/`rule_value` engine. Private equation bindings prevent
  user redefinitions of atan, log1p or the equation's letters from changing mathematical constants. No second
  transcendental series evaluator is introduced and `use prelude` is not required for definite integrals.
- `pi=4*atan(1)` is evaluated using `atan(1)=2*atan(1/3)+atan(1/7)`. Logs reduce their positive argument by powers
  of two, using `log(2)=log1p(1/2)-log1p(-1/4)`. Atan uses oddness, reciprocal arguments and subtraction of
  atan(1/2), keeping series arguments bounded away from the convergence boundary. All operations use balls.
- Atan at 0, +/-1, +/-sqrt(3), +/-1/sqrt(3) reduces exactly to rational multiples of pi. Other values remain
  named atan constants. `log(1)=0`; general logarithmic identities are not a symbolic simplification engine.
- Named-number identity includes its value callback as well as its display text. A moved user-defined rule
  named atan must not supply the value of a circular-area constant with the same printed spelling.
- Tests compare the five supplied definite values and additional range/domain cases with independent bc -l
  values at 60 places. Regression counts and the full acceptance table are in SLICE9_VERIFICATION.md.
