/* Areas under the square root of a quantity in the flowing letter.
 *
 * Newton puts a letter for the root ("pro singulis pono totidem literas", Methodus, Problem 1): here s, with
 * s^2 = q. The area under sqrt(aa + xx) is the hyperbola's, under sqrt(aa - xx) the circle's (De Analysi,
 * NATP00204): so a root of a quadratic leaves, besides an algebraic part, the area of a conic, a logarithm for
 * the hyperbola and an arc (asin) for the circle. The algebraic part is found as he finds curves whose areas are
 * known (De Analysi: "supponas sqrt(aa + xx) = z, ex calculo invenies x/sqrt(aa + xx) = y"): suppose the area
 * Q(x) s + lambda * (conic area), take its moment, and compare the terms. A root of a linear quantity is put as
 * the new flowing letter, and the area is that of a rational curve (Slice 9). Every result is put back. */
#include "nm.h"

#include <stdio.h>
#include <string.h>

static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }
static C N(int64_t v) { return c_const(qi(v)); }

static int deg_in(C p, int v) {
    int d = 0;
    for (int t = 0; t < p.nt; t++) {
        Q e = p.t[t].e[v];
        if (!q_is_int(e) || q_sign(e) < 0) nm_fail("internal: not a polynomial in the root letter");
        int64_t k; z_fits_i64(e.num, &k);
        if (k > d) d = (int)k;
    }
    return d;
}

static C coeff(C p, int v, int k) {                 /* the coefficient of v^k */
    C r = c_zero();
    for (int t = 0; t < p.nt; t++)
        if (q_cmp(p.t[t].e[v], qi(k)) == 0) { C m; m.nt = 1; m.t = arena_alloc(sizeof(CT)); m.t[0] = p.t[t]; m.t[0].e[v] = qi(0); r = c_add(r, m); }
    return r;
}

static C mono(int v, int k) { return c_pow_int(c_letter(v), k); }

/* s^2 -> q: a + b s */
static void reduce(C p, int s, C q, C *a, C *b) {
    *a = c_zero(); *b = c_zero();
    int d = deg_in(p, s);
    for (int k = 0; k <= d; k++) {
        C c = c_mul(coeff(p, s, k), c_pow_int(q, k / 2));
        if (k & 1) *b = c_add(*b, c); else *a = c_add(*a, c);
    }
}

static C diff(C p, int v) { return integ_diff_poly(p, v); }

/* a/b with b = b0 + b1 s rationalised: the parts free of s, A + B s */
static void split(R f, int x, int s, C q, R *A, R *B) {
    C n0, n1, d0, d1;
    reduce(f.num, s, q, &n0, &n1);
    reduce(f.den, s, q, &d0, &d1);
    C den = c_sub(c_mul(d0, d0), c_mul(c_mul(d1, d1), q));         /* (d0 + d1 s)(d0 - d1 s) */
    if (c_is_zero(den)) nm_fail("division by zero under the root");
    C a = c_sub(c_mul(n0, d0), c_mul(c_mul(n1, d1), q)), b = c_sub(c_mul(n1, d0), c_mul(n0, d1));
    *A = r_make(a, den); *B = r_make(b, den);
    (void)x;
}

/* the parts of f free of s: f = A + B s */
void sqrt_integral_split(R f, int s, C q, R *A, R *B) { split(f, 0, s, q, A, B); }

static int is_number(C c) { Q k; return c_is_zero(c) || c_const_value(c, &k); }

static void add_term(Integral *I, R coef, C poly, int kind) {
    AreaTerm *t = arena_alloc((size_t)(I->n + 1) * sizeof(AreaTerm));
    if (I->n) memcpy(t, I->term, (size_t)I->n * sizeof(AreaTerm));
    t[I->n].coef = coef; t[I->n].poly = poly; t[I->n].kind = kind;
    I->term = t; I->n++;
}

static void merge(Integral *I, Integral J) {
    I->rational = r_add(I->rational, J.rational);
    for (int i = 0; i < J.n; i++) add_term(I, J.term[i].coef, J.term[i].poly, J.term[i].kind);
}

/* p(x) at x = X, a quotient */
static R at(C p, int x, R X) {
    int d = deg_in(p, x);
    R r = r_from_c(c_zero());
    for (int k = d; k >= 0; k--) r = r_add(r_mul(r, X), r_from_c(coeff(p, x, k)));
    return r;
}

static int rational_sqrt(Q v, Q *r) {             /* v a square of a rational? */
    if (q_sign(v) < 0) return 0;
    Z a = z_iroot(v.num, 2), b = z_iroot(v.den, 2);
    if (z_cmp(z_mul(a, a), v.num) || z_cmp(z_mul(b, b), v.den)) return 0;
    *r = q_make(a, b);
    return 1;
}

/* A rational point (k, c) of the conic s^2 = q: x = 0 if q(0) is a square, a rational root of q, or a small
 * trial, as Newton tries small numbers for a root ("tento numeros"). */
static int conic_point(C q, int x, Q *k, Q *c) {
    Q al, be, ga, d;
    c_const_value(coeff(q, x, 2), &al); c_const_value(coeff(q, x, 1), &be); c_const_value(coeff(q, x, 0), &ga);
    if (rational_sqrt(ga, c)) { *k = qi(0); return 1; }
    Q disc = q_sub(q_mul(be, be), q_mul(q_mul(qi(4), al), ga));
    if (rational_sqrt(disc, &d)) { *k = q_div(q_sub(q_neg(be), d), q_mul(qi(2), al)); *c = qi(0); return 1; }
    for (int den = 1; den <= 12; den++)
        for (int num = -24; num <= 24; num++) {
            Q t = q_make(z_from_i64(num), z_from_i64(den)), v = q_add(q_mul(q_add(q_mul(al, t), be), t), ga);
            if (rational_sqrt(v, c) && q_sign(*c) > 0) { *k = t; return 1; }
        }
    return 0;
}

static void merge(Integral *I, Integral J);

/* Euler's substitution through a rational point of the conic: t = (s - c)/(x - k), so that s = c + t(x - k)
 * and x = (2ct - k t^2 - al k - be)/(al - t^2): x and s are quotients in t, and the area a rational one. */
static void euler(Integral *I, R B, int x, int s, C q) {
    Q k, c;
    if (!conic_point(q, x, &k, &c)) nm_fail("no rational point found on s^2 = %s; this area comes later", c_to_str(q));
    C al = coeff(q, x, 2), be = coeff(q, x, 1), K = c_const(k), Cc = c_const(c);
    char *sq = (char *)letter_name(s), *km = c_to_str(c_sub(c_letter(x), K)), *num, *name;
    num = q_sign(c) ? arena_alloc(strlen(sq) + 32) : sq;
    if (q_sign(c)) sprintf(num, "(%s - %s)", sq, q_to_str(c));
    name = arena_alloc(strlen(num) + strlen(km) + 8);
    if (q_sign(k)) sprintf(name, "(%s/(%s))", num, km); else sprintf(name, "(%s/%s)", num, km);
    int t = letter_shown(name);
    C T = c_letter(t), T2 = c_mul(T, T);
    R X = r_make(c_sub(c_sub(c_sub(c_mul(c_mul(N(2), Cc), T), c_mul(K, T2)), c_mul(al, K)), be), c_sub(al, T2));
    R Sv = r_add(r_from_c(Cc), r_mul(r_from_c(T), r_sub(X, r_from_c(K))));
    if (!r_equal(r_mul(Sv, Sv), at(q, x, X))) nm_fail("internal check failed: Euler's substitution (vitiose)");
    R dX = r_make(c_sub(c_mul(integ_diff_poly(X.num, t), X.den), c_mul(X.num, integ_diff_poly(X.den, t))), c_mul(X.den, X.den));
    R g = r_mul(r_mul(r_div(at(B.num, x, X), at(B.den, x, X)), Sv), dX);
    Integral G = integ_rational(g, t);               /* checked there: its moment is g */
    R Tx = r_div(r_sub(r_from_c(c_letter(s)), r_from_c(Cc)), r_from_c(c_sub(c_letter(x), K)));
    R alg = r_div(at(G.rational.num, t, Tx), at(G.rational.den, t, Tx)), A, Bs;
    split(alg, x, s, q, &A, &Bs);                     /* the algebraic part back in x and the root: -sqrt(x^2 + 1)/x */
    I->rational = r_add(I->rational, r_add(A, r_mul(Bs, r_from_c(c_letter(s)))));
    G.rational = r_from_c(c_zero());
    merge(I, G);
    I->euler = t + 1; I->ek = K; I->ec = Cc;
}

Integral sqrt_integral(R f, int x, int s, C q) {
    int dq = deg_in(q, x);
    for (int l = 0; l < letter_count(); l++)
        if (l != x && c_uses(q, l)) nm_fail("letters under the root besides %s come later", letter_name(x));
    if (dq < 1 || dq > 2) nm_fail("a root of degree %d in %s: its area is not a conic's; ask for the series with 'to %s^N'", dq, letter_name(x), letter_name(x));
    R A, B;
    split(f, x, s, q, &A, &B);
    Integral I; memset(&I, 0, sizeof I);
    I.var = x; I.rational = r_from_c(c_zero()); I.root = s + 1;
    if (!r_is_zero(A)) merge(&I, integ_rational(A, x));
    if (r_is_zero(B)) return I;
    if (dq == 1) {                                     /* s the new flowing letter: x = (s^2 - b)/a, dx = 2s/a ds */
        C a = coeff(q, x, 1), b = coeff(q, x, 0);
        C xs = c_div(c_sub(mono(s, 2), b), a);
        R g = r_make(integ_subst_poly(B.num, x, xs), integ_subst_poly(B.den, x, xs));
        g = r_mul(g, r_from_c(c_div(c_mul(N(2), mono(s, 2)), a)));
        Integral G = integ_rational(g, s);           /* checked there: its moment is g */
        merge(&I, G);
        return I;
    }
    /* quadratic: B q must be a polynomial P; area = Q(x) s + lambda * (area of dx/s) */
    R Pq = r_mul(B, r_from_c(q));
    if (!is_number(Pq.den) || c_uses(Pq.den, x)) { euler(&I, B, x, s, q); return I; }
    C P = c_div(Pq.num, Pq.den);
    C al = coeff(q, x, 2), be = coeff(q, x, 1), ga = coeff(q, x, 0);
    C disc = c_sub(c_mul(be, be), c_mul(c_mul(N(4), al), ga));
    if (c_is_zero(disc)) nm_fail("the root of a perfect square: write the absolute value");
    int dp = deg_in(P, x), nq = dp > 0 ? dp : 0;      /* unknowns: Q of degree dp - 1 (nq coefficients), lambda */
    int nu = nq + 1;
    Mat M = mat_new(dp + 1, nu), rhs = mat_new(dp + 1, 1);
    C qp = diff(q, x);
    for (int j = 0; j < nq; j++) {                     /* Q = x^j: Q' q + Q q'/2 */
        C e = c_add(c_mul(c_mul(N(j), mono(x, j > 0 ? j - 1 : 0)), j > 0 ? q : c_zero()), c_div(c_mul(mono(x, j), qp), N(2)));
        for (int k = 0; k <= dp; k++) M.a[k * nu + j] = r_from_c(coeff(e, x, k));
    }
    M.a[0 * nu + nq] = r_from_c(N(1));
    for (int k = 0; k <= dp; k++) rhs.a[k] = r_from_c(coeff(P, x, k));
    Mat ns, sol = mat_solve(M, rhs, &ns);
    C Qp = c_zero();
    for (int j = 0; j < nq; j++) Qp = c_add(Qp, c_mul(c_div(sol.a[j].num, sol.a[j].den), mono(x, j)));
    C lam = c_div(sol.a[nq].num, sol.a[nq].den);
    C back = c_add(c_add(c_mul(diff(Qp, x), q), c_div(c_mul(Qp, qp), N(2))), lam);
    if (!c_equal(back, P)) nm_fail("internal check failed: the assumed area does not give the curve back (vitiose)");
    I.rational = r_add(I.rational, r_from_c(c_mul(Qp, c_letter(s))));
    if (c_is_zero(lam)) return I;
    Q alq, dq2;
    c_const_value(al, &alq); c_const_value(disc, &dq2);
    if (q_sign(alq) > 0) {                           /* the hyperbola: log|2 al x + be + 2 sqrt(al) s| / sqrt(al) */
        C ra = c_radical_q(alq, 2);
        C p = c_add(c_add(c_mul(c_mul(N(2), al), c_letter(x)), be), c_mul(c_mul(N(2), ra), c_letter(s)));
        C sc = c_inv(c_mul(N(2), al));               /* a constant factor only shifts the area: x + sqrt(x^2 + 1) */
        p = c_mul(p, sc);
        /* put back: p' s = sqrt(al) p, with s' = q'/(2s) */
        C pps = c_mul(c_add(c_mul(c_mul(N(2), al), c_letter(s)), c_mul(ra, qp)), sc);
        if (!c_equal(pps, c_mul(ra, p))) nm_fail("internal check failed: the hyperbola's area (vitiose)");
        add_term(&I, r_from_c(c_div(lam, ra)), p, AREA_LOG);
    } else {                                          /* the circle: asin((-2 al x - be)/sqrt(disc)) / sqrt(-al) */
        if (q_sign(dq2) <= 0) nm_fail("the quantity under the root is negative for every %s", letter_name(x));
        C rd = c_radical_q(dq2, 2), rm = c_radical_q(q_neg(alq), 2);
        C u = c_div(c_sub(c_mul(c_mul(N(-2), al), c_letter(x)), be), rd);
        /* put back: 1 - u^2 = k q with k = -4 al/disc, and u'/sqrt(k) = sqrt(-al) */
        C k = c_div(c_mul(N(-4), al), disc);
        if (!c_equal(c_sub(N(1), c_mul(u, u)), c_mul(k, q))) nm_fail("internal check failed: the circle's area (vitiose)");
        Q kq; c_const_value(k, &kq);
        if (!c_equal(c_div(diff(u, x), c_radical_q(kq, 2)), rm)) nm_fail("internal check failed: the circle's arc (vitiose)");
        add_term(&I, r_from_c(c_div(lam, rm)), u, AREA_ASIN);
    }
    return I;
}

/* asin(u) for an exact number u: the arc of the circle, as atan(u/sqrt(1 - u^2)) or a part of pi */
C asin_constant(C u, C (*conic)(int, C)) {
    Q u2;
    if (!c_const_value(c_mul(u, u), &u2)) nm_fail("asin here takes a number whose square is rational");
    int c = q_cmp_one(u2);
    if (c > 0) nm_fail("asin of a number beyond 1: outside the circle");
    if (c == 0) {                                     /* +-pi/2: twice the arc of atan(1) */
        Q uq;
        if (!c_const_value(u, &uq)) nm_fail("internal: asin of a unit that is not a number");
        C r = c_mul(conic(AREA_ATAN, N(1)), N(2));
        return q_sign(uq) < 0 ? c_neg(r) : r;
    }
    if (c_is_zero(u)) return c_zero();
    C w = c_radical_q(q_sub(qi(1), u2), 2);
    return conic(AREA_ATAN, c_div(u, w));
}

/* F(x) at a number: s = sqrt(q(x)) first, then x */
C sqrt_integral_value(Integral a, C q, C x, C (*conic)(int, C)) {
    if (!a.root) return integ_value(a, x, conic);
    int s = a.root - 1;
    C qx = integ_subst_poly(q, a.var, x);
    Q qv;
    if (!c_const_value(qx, &qv)) nm_fail("the root at an endpoint must be of a rational number");
    if (q_sign(qv) < 0) nm_fail("the quantity under the root is negative at an endpoint");
    C sv = c_radical_q(qv, 2);
    if (a.euler) {                                    /* t at x: (s - c)/(x - k), or its limit q'(k)/(2c) at x = k */
        int t = a.euler - 1;
        C tv;
        if (c_equal(x, a.ek)) {
            if (c_is_zero(a.ec)) nm_fail("an end at the point of Euler's substitution, on the conic's root: later");
            tv = c_div(integ_subst_poly(integ_diff_poly(q, a.var), a.var, x), c_mul(N(2), a.ec));
        } else tv = c_div(c_sub(sv, a.ec), c_sub(x, a.ek));
        a.rational.num = integ_subst_poly(a.rational.num, t, tv);
        a.rational.den = integ_subst_poly(a.rational.den, t, tv);
        AreaTerm *tt = arena_alloc((size_t)(a.n + 1) * sizeof(AreaTerm));
        for (int i = 0; i < a.n; i++) { tt[i] = a.term[i]; tt[i].poly = integ_subst_poly(a.term[i].poly, t, tv); }
        a.term = tt;
    }
    Integral b = a;
    b.rational.num = integ_subst_poly(a.rational.num, s, sv);
    b.rational.den = integ_subst_poly(a.rational.den, s, sv);
    b.term = arena_alloc((size_t)(a.n + 1) * sizeof(AreaTerm));
    C y = c_zero();
    for (int i = 0; i < a.n; i++) {
        AreaTerm t = a.term[i];
        C arg = integ_subst_poly(integ_subst_poly(t.poly, s, sv), a.var, x);
        C v = t.kind == AREA_ASIN ? asin_constant(arg, conic) : conic(t.kind, arg);
        y = c_add(y, c_mul(c_div(t.coef.num, t.coef.den), v));
    }
    b.n = 0;
    return c_add(y, integ_value(b, x, conic));
}
