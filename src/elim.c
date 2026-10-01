/* Systems of equations, by extermination (Newton's "extermino literas"): the resultant of two equations in a
 * letter is an equation without it; letters are removed one by one until one remains. Its roots are found
 * (rational ones exactly, real ones as surds with certified brackets, complex ones of quadratic factors with i),
 * carried back through the other equations, and every solution is substituted into the whole system before it
 * is accepted, since extermination can bring in roots that belong to no solution. */
#include "nm.h"

#include <stdio.h>
#include <string.h>

static Q q0(void) { return q_from_z(z_zero()); }
static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }

/* ---------------- polynomials in one letter, coefficients quantities in the others ---------------- */

UP up_alloc(int deg) {
    UP p; p.deg = deg; p.c = arena_alloc((size_t)(deg + 2) * sizeof(C));
    for (int i = 0; i <= deg; i++) p.c[i] = c_zero();
    return p;
}
UP up_trim(UP p) { while (p.deg >= 0 && c_is_zero(p.c[p.deg])) p.deg--; return p; }

static int64_t whole_exp(Q e) {
    int64_t v;
    if (!q_is_int(e) || !z_fits_i64(e.num, &v) || v < 0 || v > 10000) nm_fail("an equation must have whole positive powers of its unknowns");
    return v;
}

UP up_from(C a, int v) {
    int deg = 0;
    for (int t = 0; t < a.nt; t++) { int64_t e = whole_exp(a.t[t].e[v]); if (e > deg) deg = (int)e; }
    UP p = up_alloc(deg);
    for (int t = 0; t < a.nt; t++) {
        int64_t e = whole_exp(a.t[t].e[v]);
        C one; one.nt = 1; one.t = arena_alloc(sizeof(CT)); one.t[0] = a.t[t]; one.t[0].e[v] = q0();
        p.c[e] = c_add(p.c[e], one);
    }
    return up_trim(p);
}

UP up_sub(UP a, UP b) {
    UP r = up_alloc(a.deg > b.deg ? a.deg : b.deg);
    for (int i = 0; i <= r.deg; i++) r.c[i] = c_sub(i <= a.deg ? a.c[i] : c_zero(), i <= b.deg ? b.c[i] : c_zero());
    return up_trim(r);
}

/* remainder of a by b over the field of numbers and surds (b's leading coefficient must be invertible) */
static UP up_rem(UP a, UP b) {
    C inv = c_inv(b.c[b.deg]);
    while (a.deg >= b.deg && a.deg >= 0) {
        int sh = a.deg - b.deg;
        C f = c_mul(a.c[a.deg], inv);
        UP m = up_alloc(b.deg + sh);
        for (int i = 0; i <= b.deg; i++) m.c[i + sh] = c_mul(b.c[i], f);
        UP na = up_sub(a, m);
        if (na.deg >= a.deg) { na.c[a.deg] = c_zero(); na = up_trim(na); }
        a = na;
    }
    return a;
}

static UP up_gcd(UP a, UP b) {
    while (b.deg >= 0) { UP r = up_rem(a, b); a = b; b = r; }
    return a;
}

/* ---------------- the resultant: a determinant without division (Berkowitz) ---------------- */

/* det of an n x n matrix of quantities, by Berkowitz's recurrence: only sums and products, so letters and surds
 * in the entries cause no trouble */
C det_nodiv(C *m, int n) {
    if (n == 0) return c_const(qi(1));
    C *V = arena_alloc((size_t)(n + 2) * sizeof(C));
    V[0] = c_const(qi(1)); V[1] = c_neg(m[0]);
    int len = 2;
    for (int r = 1; r < n; r++) {
        C *q = arena_alloc((size_t)(r + 2) * sizeof(C));
        q[0] = c_const(qi(1));
        q[1] = c_neg(m[r * n + r]);
        C *w = arena_alloc((size_t)r * sizeof(C)), *w2 = arena_alloc((size_t)r * sizeof(C));
        for (int i = 0; i < r; i++) w[i] = m[i * n + r];                 /* the column above the corner */
        for (int k = 0; k < r; k++) {
            C s = c_zero();
            for (int i = 0; i < r; i++) s = c_add(s, c_mul(m[r * n + i], w[i]));  /* the row left of the corner */
            q[k + 2] = c_neg(s);
            for (int i = 0; i < r; i++) {
                C t = c_zero();
                for (int j = 0; j < r; j++) t = c_add(t, c_mul(m[i * n + j], w[j]));
                w2[i] = t;
            }
            C *tmp = w; w = w2; w2 = tmp;
        }
        C *NV = arena_alloc((size_t)(r + 3) * sizeof(C));
        for (int i = 0; i <= r + 1; i++) {
            C s = c_zero();
            for (int j = 0; j <= i && j < len; j++) s = c_add(s, c_mul(q[i - j], V[j]));
            NV[i] = s;
        }
        V = NV; len = r + 2;
    }
    return n % 2 ? c_neg(V[n]) : V[n];
}

C elim_resultant(C A, C B, int v) {
    UP a = up_from(A, v), b = up_from(B, v);
    if (a.deg < 0 || b.deg < 0) return c_zero();
    if (a.deg == 0) return c_pow_int(a.c[0], b.deg);
    if (b.deg == 0) return c_pow_int(b.c[0], a.deg);
    int n = a.deg + b.deg;
    C *m = arena_alloc((size_t)n * n * sizeof(C));
    for (int i = 0; i < n * n; i++) m[i] = c_zero();
    for (int i = 0; i < b.deg; i++) for (int j = 0; j <= a.deg; j++) m[i * n + i + j] = a.c[a.deg - j];
    for (int i = 0; i < a.deg; i++) for (int j = 0; j <= b.deg; j++) m[(b.deg + i) * n + i + j] = b.c[b.deg - j];
    return det_nodiv(m, n);
}

/* ---------------- roots of one equation with rational coefficients ---------------- */

typedef struct { int deg; Q *c; } RQ;          /* rational coefficients */

static RQ rq_trim(RQ p) { while (p.deg >= 0 && q_sign(p.c[p.deg]) == 0) p.deg--; return p; }
static RQ rq_alloc(int deg) { RQ p; p.deg = deg; p.c = arena_alloc((size_t)(deg + 2) * sizeof(Q)); for (int i = 0; i <= deg; i++) p.c[i] = q0(); return p; }
static Q rq_at(RQ p, Q x) { Q a = q0(); for (int i = p.deg; i >= 0; i--) a = q_add(q_mul(a, x), p.c[i]); return a; }
static RQ rq_deriv(RQ p) { RQ d = rq_alloc(p.deg - 1); for (int i = 1; i <= p.deg; i++) d.c[i - 1] = q_mul(p.c[i], qi(i)); return rq_trim(d); }
static void rq_divmod(RQ a, RQ b, RQ *q, RQ *r) {
    RQ qq = rq_alloc(a.deg - b.deg >= 0 ? a.deg - b.deg : 0);
    a = rq_trim(a);
    while (a.deg >= b.deg && a.deg >= 0) {
        int sh = a.deg - b.deg;
        Q f = q_div(a.c[a.deg], b.c[b.deg]);
        qq.c[sh] = q_add(qq.c[sh], f);
        RQ na = rq_alloc(a.deg);
        for (int i = 0; i <= a.deg; i++) na.c[i] = a.c[i];
        for (int i = 0; i <= b.deg; i++) na.c[i + sh] = q_sub(na.c[i + sh], q_mul(b.c[i], f));
        na.c[a.deg] = q0();
        a = rq_trim(na);
    }
    *q = rq_trim(qq); *r = a;
}
static RQ rq_gcd(RQ a, RQ b) { while (b.deg >= 0) { RQ q, r; rq_divmod(a, b, &q, &r); a = b; b = r; } return a; }

/* the number of sign changes of a Sturm sequence at x */
static int sturm_changes(RQ *s, int ns, Q x) {
    int ch = 0, last = 0;
    for (int i = 0; i < ns; i++) {
        int sg = q_sign(rq_at(s[i], x));
        if (sg == 0) continue;
        if (last && sg != last) ch++;
        last = sg;
    }
    return ch;
}

static int root_counter;

/* a real root isolated in (lo, hi): narrowed by halving, then kept as a surd with its certified bracket */
static C real_root(RQ p, Q lo, Q hi, char *name_out) {
    for (int it = 0; it < 60; it++) {
        Q mid = q_div(q_add(lo, hi), qi(2));
        Q fm = rq_at(p, mid);
        if (q_sign(fm) == 0) return c_const(mid);
        if (q_sign(rq_at(p, lo)) * q_sign(fm) < 0) hi = mid; else lo = mid;
    }
    Poly pp; pp.deg = p.deg; pp.var = "y"; pp.c = p.c;
    Q mid = q_div(q_add(lo, hi), qi(2));
    Root *r = root_new(pp, mid, 20);
    root_refine(r, 20);
    if (!r->certified) nm_fail("internal: an isolated root could not be certified");
    char nm[32];
    snprintf(nm, sizeof nm, "r_%d", ++root_counter);
    if (name_out) strcpy(name_out, nm);
    return c_surd_from_root(nm, r);
}

static RQ rq_rational_roots(RQ p, C *out, int *nout, int max) {
    /* rational roots, p/q with p | a0 and q | an, for coefficients of moderate size */
    Q l = qi(1);
    for (int i = 0; i <= p.deg; i++) if (q_sign(p.c[i])) { Z gg = z_gcd(l.num, p.c[i].den), qq, rr; z_divmod(z_mul(l.num, p.c[i].den), gg, &qq, &rr); l = q_from_z(qq); }
    int64_t a0 = 0, an = 0;
    int zero_root = q_sign(p.c[0]) == 0;
    if (zero_root) {
        if (*nout < max) out[(*nout)++] = c_zero();
        RQ q, r; RQ x = rq_alloc(1); x.c[1] = qi(1); rq_divmod(p, x, &q, &r); p = q;
    }
    z_fits_i64(z_abs(q_mul(p.c[0], l).num), &a0); z_fits_i64(z_abs(q_mul(p.c[p.deg], l).num), &an);
    if (a0 > 0 && an > 0 && a0 < 1000000000000LL && an < 1000000000000LL && p.deg > 0) {
        for (int64_t pd = 1; pd * pd <= a0 && p.deg > 0; pd++) {
            if (a0 % pd) continue;
            int64_t pds[2] = {pd, a0 / pd};
            for (int u = 0; u < 2; u++)
                for (int64_t qd = 1; qd * qd <= an && p.deg > 0; qd++) {
                    if (an % qd) continue;
                    int64_t qds[2] = {qd, an / qd};
                    for (int w = 0; w < 2; w++)
                        for (int sg = -1; sg <= 1; sg += 2) {
                            if (p.deg <= 0) break;
                            Q v = q_make(z_from_i64(sg * pds[u]), z_from_i64(qds[w]));
                            if (q_sign(rq_at(p, v)) == 0) {
                                if (*nout < max) out[(*nout)++] = c_const(v);
                                RQ lin = rq_alloc(1); lin.c[1] = qi(1); lin.c[0] = q_neg(v);
                                RQ q, r; rq_divmod(p, lin, &q, &r); p = q;
                            }
                        }
                }
        }
    }
    return p;
}

/* the roots of p (rational coefficients): rational ones exactly, other real ones as surds, and the complex pair of
 * a remaining quadratic with i. Returns how many roots were left unresolved (complex roots of higher factors). */
static int rq_roots(RQ p, C *out, int *nout, int max) {
    p = rq_trim(p);
    if (p.deg <= 0) return 0;
    RQ g = rq_gcd(p, rq_deriv(p));                 /* the square-free part: each root once */
    if (g.deg > 0) { RQ q, r; rq_divmod(p, g, &q, &r); p = q; }
    p = rq_rational_roots(p, out, nout, max);
    if (p.deg <= 0) return 0;
    if (p.deg == 2) {                              /* a quadratic: the formula, with i when needed */
        Q A = p.c[2], B = p.c[1], Cc = p.c[0];
        Q disc = q_sub(q_mul(B, B), q_mul(qi(4), q_mul(A, Cc)));
        C s = c_radical_q(disc, 2);
        C mB = c_const(q_neg(B)), twoA = c_const(q_mul(qi(2), A));
        if (*nout < max) out[(*nout)++] = c_div(c_add(mB, s), twoA);
        if (*nout < max) out[(*nout)++] = c_div(c_sub(mB, s), twoA);
        return 0;
    }
    /* real roots by Sturm's sequence and halving */
    RQ s[256]; int ns = 0;
    s[ns++] = p; s[ns++] = rq_deriv(p);
    while (ns < 256 && s[ns - 1].deg > 0) {
        RQ q, r; rq_divmod(s[ns - 2], s[ns - 1], &q, &r);
        if (r.deg < 0) break;
        for (int i = 0; i <= r.deg; i++) r.c[i] = q_neg(r.c[i]);
        s[ns++] = r;
    }
    Q bound = qi(1);                               /* every root lies within 1 + max |a_i / a_n| */
    for (int i = 0; i < p.deg; i++) { Q v = q_div(p.c[i], p.c[p.deg]); if (q_sign(v) < 0) v = q_neg(v); bound = q_add(bound, v); }
    Q stack_lo[512], stack_hi[512]; int sp = 0, found = 0;
    stack_lo[sp] = q_neg(bound); stack_hi[sp] = bound; sp++;
    while (sp > 0) {
        sp--;
        Q lo = stack_lo[sp], hi = stack_hi[sp];
        int k = sturm_changes(s, ns, lo) - sturm_changes(s, ns, hi);
        if (k == 0) continue;
        if (k == 1) { if (*nout < max) out[(*nout)++] = real_root(p, lo, hi, NULL); found++; continue; }
        Q mid = q_div(q_add(lo, hi), qi(2));
        if (q_sign(rq_at(p, mid)) == 0) { if (*nout < max) out[(*nout)++] = c_const(mid); found++; mid = q_add(mid, q_div(q_sub(hi, lo), qi(1000))); }
        if (sp + 2 >= 512) nm_fail("too many roots to separate");
        stack_lo[sp] = lo; stack_hi[sp] = mid; sp++;
        stack_lo[sp] = mid; stack_hi[sp] = hi; sp++;
    }
    return p.deg - found;                          /* complex roots of a factor above degree two */
}

/* roots of a polynomial whose coefficients may hold surds: degree one and two exactly */
static int up_roots(UP p, C *out, int *nout, int max) {
    p = up_trim(p);
    if (p.deg <= 0) return 0;
    int rational = 1;
    for (int i = 0; i <= p.deg; i++) { Q k; if (!c_const_value(p.c[i], &k)) rational = 0; }
    if (rational) {
        RQ r = rq_alloc(p.deg);
        for (int i = 0; i <= p.deg; i++) c_const_value(p.c[i], &r.c[i]);
        return rq_roots(r, out, nout, max);
    }
    if (c_is_zero(p.c[0])) {                       /* the root 0, then the rest */
        if (*nout < max) out[(*nout)++] = c_zero();
        UP q = up_alloc(p.deg - 1);
        for (int i = 1; i <= p.deg; i++) q.c[i - 1] = p.c[i];
        return up_roots(up_trim(q), out, nout, max);
    }
    if (p.deg == 1) { if (*nout < max) out[(*nout)++] = c_neg(c_div(p.c[0], p.c[1])); return 0; }
    if (p.deg == 2) {
        C A = p.c[2], B = p.c[1], Cc = p.c[0];
        C disc = c_sub(c_mul(B, B), c_scale(c_mul(A, Cc), qi(4)));
        C s;
        Q k;
        if (c_const_value(disc, &k)) s = c_radical_q(k, 2);
        else if (!c_pow_q(disc, q_make(z_from_i64(1), z_from_i64(2)), &s)) {
            if (c_has_plain(disc)) return 2;          /* sqrt of a sum of letters: later */
            s = c_radical_c(disc, 2);
        }
        C twoA = c_scale(A, qi(2));
        if (*nout < max) out[(*nout)++] = c_div(c_sub(s, B), twoA);
        if (*nout < max) out[(*nout)++] = c_div(c_sub(c_neg(B), s), twoA);
        return 0;
    }
    return p.deg;
}

/* ---------------- solving a system ---------------- */

static C subst(C a, int v, C val) {               /* a with the letter v replaced by val */
    UP p = up_from(a, v);
    C out = c_zero();
    for (int i = p.deg; i >= 0; i--) out = c_add(c_mul(out, val), p.c[i]);
    return out;
}

static int uses_any(C a, int *vars, int nv) {
    for (int i = 0; i < nv; i++) if (c_uses(a, vars[i])) return 1;
    return 0;
}

/* solutions of eqs = 0 in the letters vars[0..nv); each solution is an array of nv values in that order */
static void solve_rec(C *eqs, int ne, int *vars, int nv, Solutions *S) {
    /* drop the vanishing equations; a nonzero one without unknowns means no solution here */
    C *E = arena_alloc((size_t)(ne + 1) * sizeof(C));
    int k = 0;
    for (int i = 0; i < ne; i++) {
        if (c_is_zero(eqs[i])) continue;
        if (!uses_any(eqs[i], vars, nv)) {
            if (c_has_plain(eqs[i]))
                nm_fail("an equation left without unknowns still holds letters (%s); solutions for all their values come later", c_to_str(eqs[i]));
            return;
        }
        E[k++] = eqs[i];
    }
    if (k == 0) { S->curve = 1; return; }          /* nothing left to fix the unknowns */
    int v = vars[nv - 1];
    int piv = -1, pdeg = 0;
    for (int i = 0; i < k; i++) {
        if (!c_uses(E[i], v)) continue;
        int d = up_from(E[i], v).deg;
        if (piv < 0 || d < pdeg) { piv = i; pdeg = d; }
    }
    if (piv < 0) { S->curve = 1; return; }          /* v is free: the solutions are not isolated */
    if (nv == 1) {                                  /* one letter: the roots of what the equations share */
        UP g = up_from(E[piv], v);
        for (int i = 0; i < k; i++) if (i != piv) g = up_gcd(g, up_from(E[i], v));
        C roots[256]; int nroots = 0;
        S->unresolved += up_roots(g, roots, &nroots, 256);
        for (int i = 0; i < nroots && S->n < NM_MAXSOL; i++) {
            S->val[S->n] = arena_alloc(sizeof(C));
            S->val[S->n][0] = roots[i];
            S->n++;
        }
        return;
    }
    /* exterminate v: the resultant of the pivot with every other equation that holds v */
    C *R = arena_alloc((size_t)(k + 1) * sizeof(C));
    int nr = 0;
    for (int i = 0; i < k; i++) {
        if (i == piv) continue;
        R[nr++] = c_uses(E[i], v) ? elim_resultant(E[piv], E[i], v) : E[i];
    }
    if (nr == 0) { S->curve = 1; return; }
    Solutions sub;
    memset(&sub, 0, sizeof sub);
    sub.nv = nv - 1;
    sub.val = arena_alloc(NM_MAXSOL * sizeof(C *));
    solve_rec(R, nr, vars, nv - 1, &sub);
    S->curve |= sub.curve;
    S->unresolved += sub.unresolved;
    /* carry each partial solution back: the equations with it put in share a factor in v */
    for (int s2 = 0; s2 < sub.n; s2++) {
        UP g; g.deg = -1; g.c = NULL;
        for (int i = 0; i < k; i++) {
            C e = E[i];
            for (int j = 0; j < nv - 1; j++) e = subst(e, vars[j], sub.val[s2][j]);
            if (c_is_zero(e)) continue;
            UP u = up_from(e, v);
            g = g.deg < 0 ? u : up_gcd(g, u);
            if (g.deg == 0) break;
        }
        if (g.deg == 0) continue;                   /* brought in by the extermination: no v fits */
        if (g.deg < 0) { S->curve = 1; continue; }  /* every v fits */
        C roots[64]; int nroots = 0;
        S->unresolved += up_roots(g, roots, &nroots, 64);
        for (int i = 0; i < nroots && S->n < NM_MAXSOL; i++) {
            S->val[S->n] = arena_alloc((size_t)nv * sizeof(C));
            for (int j = 0; j < nv - 1; j++) S->val[S->n][j] = sub.val[s2][j];
            S->val[S->n][nv - 1] = roots[i];
            S->n++;
        }
    }
}

/* Solve eqs[0..ne) = 0 for the letters vars[0..nv). Every solution returned has been substituted back. */
Solutions elim_solve(C *eqs, int ne, int *vars, int nv) {
    if (nv < 1 || nv > 6) nm_fail("between 1 and 6 unknowns in this version");
    Solutions raw;
    memset(&raw, 0, sizeof raw);
    raw.nv = nv;
    raw.val = arena_alloc(NM_MAXSOL * sizeof(C *));
    solve_rec(eqs, ne, vars, nv, &raw);
    Solutions S = raw;
    S.n = 0;
    S.val = arena_alloc(NM_MAXSOL * sizeof(C *));
    for (int s = 0; s < raw.n; s++) {               /* the check: every equation must vanish */
        int ok = 1;
        for (int i = 0; i < ne && ok; i++) {
            C e = eqs[i];
            for (int j = 0; j < nv; j++) e = subst(e, vars[j], raw.val[s][j]);
            if (!c_is_zero(e)) ok = 0;
        }
        if (!ok) { S.rejected++; continue; }
        int dup = 0;
        for (int t = 0; t < S.n && !dup; t++) {
            int same = 1;
            for (int j = 0; j < nv; j++) if (!c_equal(S.val[t][j], raw.val[s][j])) same = 0;
            dup = same;
        }
        if (!dup) S.val[S.n++] = raw.val[s];
    }
    return S;
}

/* the roots of one polynomial in the letter v (numbers, surds and letters in its coefficients) */
int elim_roots(C p, int v, C *out, int max, int *unresolved) {
    int n = 0;
    *unresolved = up_roots(up_from(p, v), out, &n, max);
    return n;
}


/* Rational linear factors followed by exact quadratic factors. Kronecker's
 * three integer samples find quadratic factors missed by the root theorem.
 * Every candidate is divided back; bounded searches never assert irreducibility. */
static RQ rq_from_c(C p, int v) {
    UP a = up_from(p,v); RQ r = rq_alloc(a.deg);
    for (int i=0;i<=a.deg;i++)
        if (!c_const_value(a.c[i], &r.c[i])) nm_fail("letters in a nonlinear conic factor: sign unknown; later");
    return r;
}
static C rq_to_c(RQ p, int v) {
    C a=c_zero(), x=c_letter(v);
    for (int i=p.deg;i>=0;i--) a=c_add(c_mul(a,x),c_const(p.c[i]));
    return a;
}
static void small_factors(RQ p, int v, C *out, int *n, int max) {
    if (p.deg<=0) return;
    if (*n>=max) nm_fail("too many conic factors");
    if (p.deg==2) {
        Q disc=q_sub(q_mul(p.c[1],p.c[1]),q_mul(qi(4),q_mul(p.c[2],p.c[0]))),root;
        if(q_sign(disc)>=0 && q_root_exact(disc,2,&root)) {
            if(*n+2>max) nm_fail("too many conic factors");
            out[(*n)++]=c_sub(c_letter(v),c_const(q_div(q_sub(q_neg(p.c[1]),root),q_mul(qi(2),p.c[2]))));
            out[(*n)++]=c_sub(c_letter(v),c_const(q_div(q_add(q_neg(p.c[1]),root),q_mul(qi(2),p.c[2]))));
            return;
        }
    }
    if (p.deg<=2) { out[(*n)++]=rq_to_c(p,v); return; }
    C roots[256]; int nr=0;
    RQ rem=rq_rational_roots(p,roots,&nr,256);
    for (int i=0;i<nr;i++) {
        if (*n>=max) nm_fail("too many conic factors");
        out[(*n)++]=c_sub(c_letter(v),roots[i]);
    }
    if (nr) { small_factors(rem,v,out,n,max); return; }
    if (p.deg>=4) {
        Z l=z_from_i64(1);
        for (int i=0;i<=p.deg;i++) {
            Z g=z_gcd(l,p.c[i].den),q,r;
            z_divmod(z_mul(l,p.c[i].den),g,&q,&r); l=q;
        }
        RQ ip=rq_alloc(p.deg);
        Z content=z_zero();
        for(int i=0;i<=p.deg;i++) { ip.c[i]=q_mul(p.c[i],q_from_z(l)); content=z_gcd(content,ip.c[i].num); }
        for(int i=0;i<=p.deg;i++) ip.c[i]=q_div(ip.c[i],q_from_z(content));
        Z ds[3][256]; int nd[3];
        for(int i=0;i<3;i++) nd[i]=z_divisors(z_abs(rq_at(ip,qi(i-1)).num),ds[i],256);
        int tried=0;
        for(int i=0;i<nd[0];i++) for(int j=0;j<nd[1];j++) for(int k=0;k<nd[2];k++)
            for(int s=-1;s<=1;s+=2) for(int t=-1;t<=1;t+=2) {
                if(++tried>100000) nm_fail("quadratic factor search bound reached; later: Rothstein-Trager");
                Q fm=q_from_z(ds[0][i]), f0=q_mul(qi(s),q_from_z(ds[1][j])), fp=q_mul(qi(t),q_from_z(ds[2][k]));
                RQ f=rq_alloc(2);
                f.c[0]=f0; f.c[1]=q_div(q_sub(fp,fm),qi(2));
                f.c[2]=q_sub(q_div(q_add(fp,fm),qi(2)),f0);
                if(!q_sign(f.c[2]) || !q_is_int(f.c[1]) || !q_is_int(f.c[2])) continue;
                RQ q,r; rq_divmod(p,f,&q,&r);
                if(r.deg<0) { out[(*n)++]=rq_to_c(f,v); small_factors(q,v,out,n,max); return; }
            }
    }
    nm_fail("remaining degree-%d factor has no supported linear/quadratic split; later: Rothstein-Trager",p.deg);
}
int elim_conic_factors(C p, int v, C *out, int max) {
    int n=0; small_factors(rq_from_c(p,v),v,out,&n,max); return n;
}
int elim_has_root_closed(C p, int v, Q lo, Q hi) {
    RQ a=rq_from_c(p,v);
    if(q_cmp(lo,hi)>0) { Q t=lo;lo=hi;hi=t; }
    if(!q_sign(rq_at(a,lo)) || !q_sign(rq_at(a,hi))) return 1;
    if(a.deg<=0) return 0;
    RQ s[256]; int n=0; s[n++]=a; s[n++]=rq_deriv(a);
    while(s[n-1].deg>0) {
        if(n>=256) nm_fail("pole check degree limit");
        RQ q,r; rq_divmod(s[n-2],s[n-1],&q,&r);
        if(r.deg<0) break;
        for(int i=0;i<=r.deg;i++) r.c[i]=q_neg(r.c[i]);
        s[n++]=r;
    }
    return sturm_changes(s,n,lo)!=sturm_changes(s,n,hi);
}
