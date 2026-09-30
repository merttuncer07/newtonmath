# newtonmath: decisions log
Each decision gets: the question, the generic answer, what Newton did (sources), a recommendation, then Mert's
decision. Nothing is built before Mert decides.

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
