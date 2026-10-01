# Brief: a programming language whose only job is computing mathematics

## What I want to build (Mert, fill in / correct)
- Purpose: a language designed from scratch for mathematical computation only (no general-purpose programming).
- Users: <who writes in it? mathematicians, engineers, AI agents?>
- Kinds of mathematics: <exact algebra, numerics, series, ODEs, proofs/certificates, ...?>
- What existing tools (Mathematica, Julia, Python+sympy, Lean, ...) do badly that this should do well: <...>
- Constraints: <interpreted/compiled, runs where, exactness vs speed>

## Questions for Newton
1. What should the primitive objects and operations be? What is the smallest set from which the rest follows
   (cf. his reduction of analysis to two problems, and to operations on series)?
2. Notation first or meaning first? How should notation be chosen so that calculation becomes mechanical?
3. How should the language treat exactness vs. approximation (he computed series to many places and checked results a second way)?
4. How should a result be checked inside the language (a second route, a certificate)?
5. What makes others adopt a notation? (His dot notation lost to Leibniz's d/dx.)

## How to use this package
- agent/: the Newton decision policy (CONSTITUTION.md, RUNTIME.md, the Newtonbot skill). Decide as it says; outputs plain, no persona.
- episodes/MATH_EPISODES.md: 50 dated episodes of his mathematical decisions (from general knowledge, flagged; not verified
  against the texts).
- texts/: seven primary texts on how he built and presented his method; open them when a question needs the source.
