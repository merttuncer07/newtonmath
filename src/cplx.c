/* Values of series at complex arguments, and of their derivatives.
 *
 * A function known by a linear fluxional equation has a term rule: D(n) c_n = sum_t N_t(n) c_{n-t} - h_{n-s}.
 * Its value, or the value of its j-th fluxion, at any point A inside the circle of convergence is the sum of the
 * terms, carried as balls, with the rest bounded from the rule. The bound needs only an upper bound r of |A|,
 * so a complex A costs nothing more. As with logarithms from 0.1 and 0.2, a point near the edge is refused with
 * the reason, not guessed at. */
#include "nm.h"

#include <string.h>

static Q q0(void) { return q_from_z(z_zero()); }
static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }
static Q qabs(Q a) { return q_sign(a) < 0 ? q_neg(a) : a; }

/* ---------------- complex balls ---------------- */

CBall cb_from_q(Q re, Q im, int64_t prec) { CBall z; z.re = b_from_q(re, prec); z.im = b_from_q(im, prec); return z; }
CBall cb_add(CBall a, CBall b, int64_t prec) { CBall z; z.re = b_add(a.re, b.re, prec); z.im = b_add(a.im, b.im, prec); return z; }
CBall cb_sub(CBall a, CBall b, int64_t prec) { CBall z; z.re = b_sub(a.re, b.re, prec); z.im = b_sub(a.im, b.im, prec); return z; }

CBall cb_mul(CBall a, CBall b, int64_t prec) {
    CBall z;
    z.re = b_sub(b_mul(a.re, b.re, prec), b_mul(a.im, b.im, prec), prec);
    z.im = b_add(b_mul(a.re, b.im, prec), b_mul(a.im, b.re, prec), prec);
    return z;
}

CBall cb_div(CBall a, CBall b, int64_t prec) {
    Ball den = b_add(b_mul(b.re, b.re, prec), b_mul(b.im, b.im, prec), prec);
    CBall z;
    z.re = b_div(b_add(b_mul(a.re, b.re, prec), b_mul(a.im, b.im, prec), prec), den, prec);
    z.im = b_div(b_sub(b_mul(a.im, b.re, prec), b_mul(a.re, b.im, prec), prec), den, prec);
    return z;
}

/* an upper bound of |value|, as a rational, rounded up to 20 significant digits: a bound only needs to be a
 * bound, and a short one keeps the fraction from being reduced at every term (the gcd of 1000-digit numbers) */
static Q ball_upper(Ball b) {
    Z top = z_add(z_abs(b.m), b.r);
    int64_t e = b.e, d = z_digits(top);
    if (d > 20) {
        Z q, r;
        z_divmod(top, z_pow10(d - 20), &q, &r);
        top = z_add(q, z_from_i64(1));
        e += d - 20;
    }
    return e >= 0 ? q_from_z(z_mul_pow10(top, e)) : q_make(top, z_pow10(-e));
}

Q cb_abs_upper(CBall a) {                         /* sqrt(re^2 + im^2), rounded up */
    Q x = ball_upper(a.re), y = ball_upper(a.im);
    Q s = q_add(q_mul(x, x), q_mul(y, y));
    if (q_sign(y) == 0) return x;
    if (q_sign(x) == 0) return y;
    /* ceil(sqrt(s) * 10^20) / 10^20 */
    Z n = z_mul_pow10(s.num, 40), q, r;
    z_divmod(n, s.den, &q, &r);
    Z root = z_add(z_iroot(q, 2), z_from_i64(1));
    return q_make(root, z_pow10(20));
}

static Ball widen(Ball b, Q err) {
    Z num = z_abs(err.num), den = err.den, q, rem;
    if (b.e < 0) num = z_mul_pow10(num, -b.e); else den = z_mul(den, z_pow10(b.e));
    z_divmod(num, den, &q, &rem);
    b.r = z_add(z_add(b.r, q), z_from_i64(rem.s ? 1 : 0));
    return b;
}

/* ---------------- the value from the term rule ---------------- */

static Q p_at_q(Poly p, Q x) { Q a = q0(); for (int i = p.deg; i >= 0; i--) a = q_add(q_mul(a, x), p.c[i]); return a; }

/* an upper bound for |P(n)/D(n)| valid for every n >= N (deg P <= deg D); negative if none yet */
static Q ratio_bound(Poly P, Poly D, int64_t N) {
    int d = D.deg;
    Q nn = qi(N), num = q0(), den = qabs(D.c[d]);
    for (int i = 0; i <= P.deg; i++) num = q_add(num, q_div(qabs(P.c[i]), q_pow(nn, d - i)));
    for (int i = 0; i < d; i++) den = q_sub(den, q_div(qabs(D.c[i]), q_pow(nn, d - i)));
    if (q_sign(den) <= 0) return qi(-1);
    return q_div(num, den);
}

/* the limit of sum_t |N_t/D| r^t: below 1, the terms end by shrinking geometrically */
int rule_converges(const TermRule *R, Q r) {
    Q g = q0();
    for (int t = 1; t <= R->T; t++) {
        if (R->N[t].deg > R->D.deg) return 0;
        if (R->N[t].deg == R->D.deg) g = q_add(g, q_mul(qabs(q_div(R->N[t].c[R->D.deg], R->D.c[R->D.deg])), q_pow(r, t)));
    }
    return q_cmp_one(g) < 0;
}

/* sum_{k >= j} c_k k!/(k-j)! A^(k-j): the j-th fluxion at A. seed: c_0 .. c_{nseed-1}; the rule gives the rest
 * from `start` on (where D(n) != 0 and h no longer enters). places: the absolute accuracy wanted. */
CBall rule_value(const TermRule *R, const Q *seed, int nseed, int64_t start, CBall A, int deriv, int64_t places, int64_t prec) {
    Q r = cb_abs_upper(A);
    {                                              /* round r up to a short fraction: bounds need no more */
        Z q, rem;
        z_divmod(z_mul_pow10(r.num, 15), r.den, &q, &rem);
        r = q_make(z_add(q, z_from_i64(1)), z_pow10(15));
    }
    Ball rb = b_from_q(r, prec);
    if (!rule_converges(R, r)) nm_fail("the terms shrink too slowly at this point to bound the rest; use a point nearer 0");
    int cap = nseed + 64;
    Q *c = arena_alloc((size_t)cap * sizeof(Q));
    for (int i = 0; i < nseed; i++) c[i] = seed[i];
    int ncoef = nseed;
    Q eps = q_make(z_from_i64(1), z_pow10(places + 10));
    CBall S = cb_from_q(q0(), q0(), prec), P = cb_from_q(qi(1), q0(), prec);
    CBall *term = arena_alloc((size_t)cap * sizeof(CBall));
    Q *d = arena_alloc((size_t)cap * sizeof(Q));
    Q *T = arena_alloc((size_t)cap * sizeof(Q));       /* T[m]: an upper bound of |d_m| r^m */
    Ball Rpow = b_from_q(qi(1), prec);
    for (int m = 0;; m++) {
        int k = m + deriv;
        while (k >= ncoef) {
            if (ncoef + 1 >= cap) {
                int ncap = cap * 2;
                Q *nc = arena_alloc((size_t)ncap * sizeof(Q)), *nd = arena_alloc((size_t)ncap * sizeof(Q));
                Q *nT = arena_alloc((size_t)ncap * sizeof(Q));
                CBall *nt = arena_alloc((size_t)ncap * sizeof(CBall));
                memcpy(nc, c, (size_t)ncoef * sizeof(Q)); memcpy(nd, d, (size_t)m * sizeof(Q)); memcpy(nt, term, (size_t)m * sizeof(CBall));
                memcpy(nT, T, (size_t)m * sizeof(Q));
                c = nc; d = nd; term = nt; T = nT; cap = ncap;
            }
            int n = ncoef;
            Q acc = q0();
            for (int t = 1; t <= R->T; t++) acc = q_add(acc, q_mul(p_at_q(R->N[t], qi(n)), c[n - t]));
            int mm = n - R->s;
            if (mm >= 0 && mm <= R->h.deg) acc = q_sub(acc, R->h.c[mm]);
            c[n] = q_div(acc, p_at_q(R->D, qi(n)));
            ncoef++;
        }
        Q f = c[k];                                    /* c_k k!/(k-j)! */
        for (int u = 1; u <= deriv; u++) f = q_mul(f, qi(m + u));
        d[m] = f;
        T[m] = ball_upper(b_mul(b_from_q(qabs(f), prec), Rpow, prec));
        Rpow = b_mul(Rpow, rb, prec);
        CBall t;
        if (q_sign(f)) { CBall fb = cb_from_q(f, q0(), prec); t = cb_mul(fb, P, prec); }
        else t = cb_from_q(q0(), q0(), prec);
        term[m] = t;
        S = cb_add(S, t, prec);
        P = cb_mul(P, A, prec);
        int64_t M = m + 1;
        if (M + deriv < start || M < R->T + 1) continue;
        /* the careful bound only once the last terms are small */
        int small = 1;
        for (int u = 1; u <= R->T && small; u++) {
            if (q_cmp(T[M - u], q_mul(eps, qi(1000))) > 0) small = 0;
        }
        if (!small) { if (m > 2000000) nm_fail("too many terms"); continue; }
        Q gamma = q0();
        int ok = 1;
        for (int u = 1; u <= R->T; u++) {
            Q b = ratio_bound(R->N[u], R->D, M + deriv);
            if (q_sign(b) < 0 || M + 1 - u <= 0) { ok = 0; break; }
            if (deriv) b = q_mul(b, q_pow(q_div(qi(M + 1), qi(M + 1 - u)), deriv));
            gamma = q_add(gamma, q_mul(b, q_pow(r, u)));
        }
        if (!ok || q_cmp_one(gamma) >= 0) continue;
        Q W = q0();
        for (int u = 1; u <= R->T; u++) if (q_cmp(T[M - u], W) > 0) W = T[M - u];
        Q tail = q_div(q_mul(q_mul(qi(R->T), gamma), W), q_sub(qi(1), gamma));
        if (q_cmp(tail, eps) <= 0) { S.re = widen(S.re, tail); S.im = widen(S.im, tail); return S; }
    }
}

/* ---------------- named values: exp(1), sin(1) ... as letters, known by their value ---------------- */
/* (the letters themselves live in coef.c; see c_named) */
