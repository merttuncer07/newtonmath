/* Factors of a polynomial in one letter over the rationals.
 *
 * Newton's question (Arithmetica Universalis, "De inventione divisorum"): find the divisors of an equation from
 * whole-number arithmetic. His substitution of 3, 2, 1, 0, -1, -2 grows without bound with the degree; the modern
 * answer to the same question is whole-number arithmetic modulo a prime:
 *   1. square-free parts (Yun: gcds with the fluxion),
 *   2. the factors modulo a small prime p (Berlekamp, deterministic),
 *   3. Hensel's lifting to p^k past Mignotte's bound on the coefficients of any factor,
 *   4. the true factors from products of lifted ones (Zassenhaus), each divided out exactly.
 * Every combination of modular factors is tried, so a factor that is left is proved irreducible; the product is
 * multiplied back at the end. */
#include "nm.h"

#include <string.h>

static Z zi(int64_t v) { return z_from_i64(v); }

/* ---------------- integer polynomials, dense ---------------- */

typedef struct { int d; Z *c; } ZP;               /* d = -1 for 0 */

static ZP zp_alloc(int d) { ZP a; a.d = d; a.c = arena_alloc((size_t)(d + 2) * sizeof(Z)); for (int i = 0; i <= d; i++) a.c[i] = z_zero(); return a; }
static ZP zp_trim(ZP a) { while (a.d >= 0 && a.c[a.d].s == 0) a.d--; return a; }
static ZP zp_copy(ZP a) { ZP b = zp_alloc(a.d); for (int i = 0; i <= a.d; i++) b.c[i] = a.c[i]; return b; }
static ZP zp_mul(ZP a, ZP b) {
    if (a.d < 0 || b.d < 0) return zp_alloc(-1);
    ZP c = zp_alloc(a.d + b.d);
    for (int i = 0; i <= a.d; i++) if (a.c[i].s) for (int j = 0; j <= b.d; j++) if (b.c[j].s) c.c[i + j] = z_add(c.c[i + j], z_mul(a.c[i], b.c[j]));
    return zp_trim(c);
}
static ZP zp_sub(ZP a, ZP b) {
    ZP c = zp_alloc(a.d > b.d ? a.d : b.d);
    for (int i = 0; i <= c.d; i++) c.c[i] = z_sub(i <= a.d ? a.c[i] : z_zero(), i <= b.d ? b.c[i] : z_zero());
    return zp_trim(c);
}
static ZP zp_deriv(ZP a) {
    if (a.d <= 0) return zp_alloc(-1);
    ZP c = zp_alloc(a.d - 1);
    for (int i = 1; i <= a.d; i++) c.c[i - 1] = z_mul(a.c[i], zi(i));
    return zp_trim(c);
}
static ZP zp_primitive(ZP a) {                     /* gcd of the coefficients 1, first coefficient > 0 */
    a = zp_trim(a);
    if (a.d < 0) return a;
    Z g = z_zero();
    for (int i = 0; i <= a.d; i++) g = z_gcd(g, a.c[i]);
    if (a.c[a.d].s < 0) g = z_neg(g);
    ZP b = zp_alloc(a.d);
    for (int i = 0; i <= a.d; i++) { Z q, r; z_divmod(a.c[i], g, &q, &r); b.c[i] = q; }
    return b;
}
/* a = b q over the integers? */
static int zp_divexact(ZP a, ZP b, ZP *q) {
    a = zp_trim(a);
    if (a.d < 0) { *q = a; return 1; }
    if (b.d > a.d) return 0;
    *q = zp_alloc(a.d - b.d);
    ZP r = zp_copy(a);
    for (int i = a.d; i >= b.d; i--) {
        if (!r.c[i].s) continue;
        Z qq, rr; z_divmod(r.c[i], b.c[b.d], &qq, &rr);
        if (rr.s) return 0;
        q->c[i - b.d] = qq;
        for (int j = 0; j <= b.d; j++) r.c[i - b.d + j] = z_sub(r.c[i - b.d + j], z_mul(qq, b.c[j]));
    }
    for (int i = 0; i < b.d; i++) if (r.c[i].s) return 0;
    *q = zp_trim(*q);
    return 1;
}
static ZP zp_prem(ZP a, ZP b) {                    /* lc(b)^(da-db+1) a mod b */
    ZP r = zp_copy(a);
    for (int i = r.d; i >= b.d; i--) {
        Z t = r.c[i];
        for (int j = 0; j <= r.d; j++) r.c[j] = z_mul(r.c[j], b.c[b.d]);
        if (t.s) for (int j = 0; j <= b.d; j++) r.c[i - b.d + j] = z_sub(r.c[i - b.d + j], z_mul(t, b.c[j]));
        r.c[i] = z_zero();
    }
    return zp_trim(r);
}
static ZP zp_gcd(ZP a, ZP b) {                     /* primitive, by the primitive remainder sequence */
    a = zp_primitive(a); b = zp_primitive(b);
    if (a.d < b.d) { ZP t = a; a = b; b = t; }
    while (b.d >= 0) {
        ZP r = zp_primitive(zp_prem(a, b));
        a = b; b = r;
    }
    return a;
}

/* ---------------- polynomials modulo a small prime ---------------- */

typedef struct { int d; int64_t *c; } MP;

static MP mp_alloc(int d) { MP a; a.d = d; a.c = arena_alloc((size_t)(d + 2) * sizeof(int64_t)); for (int i = 0; i <= d; i++) a.c[i] = 0; return a; }
static MP mp_trim(MP a) { while (a.d >= 0 && a.c[a.d] == 0) a.d--; return a; }
static MP mp_copy(MP a) { MP b = mp_alloc(a.d); for (int i = 0; i <= a.d; i++) b.c[i] = a.c[i]; return b; }
static int64_t md(int64_t a, int64_t p) { a %= p; return a < 0 ? a + p : a; }
static int64_t inv_mod(int64_t a, int64_t p) {
    int64_t t = 0, nt = 1, r = p, nr = md(a, p);
    while (nr) { int64_t q = r / nr, x; x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
    return md(t, p);
}
static int64_t zmodp(Z a, int64_t p) { int64_t v; z_fits_i64(z_mod(a, zi(p)), &v); return md(v, p); }
static MP mp_from(ZP a, int64_t p) { MP m = mp_alloc(a.d); for (int i = 0; i <= a.d; i++) m.c[i] = zmodp(a.c[i], p); return mp_trim(m); }
static MP mp_mul(MP a, MP b, int64_t p) {
    if (a.d < 0 || b.d < 0) return mp_alloc(-1);
    MP c = mp_alloc(a.d + b.d);
    for (int i = 0; i <= a.d; i++) if (a.c[i]) for (int j = 0; j <= b.d; j++) c.c[i + j] = (c.c[i + j] + a.c[i] * b.c[j]) % p;
    return mp_trim(c);
}
static MP mp_sub(MP a, MP b, int64_t p) {
    MP c = mp_alloc(a.d > b.d ? a.d : b.d);
    for (int i = 0; i <= c.d; i++) c.c[i] = md((i <= a.d ? a.c[i] : 0) - (i <= b.d ? b.c[i] : 0), p);
    return mp_trim(c);
}
static MP mp_add(MP a, MP b, int64_t p) {
    MP c = mp_alloc(a.d > b.d ? a.d : b.d);
    for (int i = 0; i <= c.d; i++) c.c[i] = md((i <= a.d ? a.c[i] : 0) + (i <= b.d ? b.c[i] : 0), p);
    return mp_trim(c);
}
static void mp_divmod(MP a, MP b, int64_t p, MP *q, MP *r) {
    a = mp_copy(a);
    int64_t il = inv_mod(b.c[b.d], p);
    *q = mp_alloc(a.d - b.d >= 0 ? a.d - b.d : 0);
    for (int i = a.d; i >= b.d; i--) {
        int64_t f = a.c[i] * il % p;
        if (!f) continue;
        q->c[i - b.d] = f;
        for (int j = 0; j <= b.d; j++) a.c[i - b.d + j] = md(a.c[i - b.d + j] - f * b.c[j], p);
    }
    *q = mp_trim(*q); *r = mp_trim(a);
}
static MP mp_monic(MP a, int64_t p) { int64_t il = inv_mod(a.c[a.d], p); MP b = mp_alloc(a.d); for (int i = 0; i <= a.d; i++) b.c[i] = a.c[i] * il % p; return b; }
static MP mp_gcd(MP a, MP b, int64_t p) {
    while (b.d >= 0) { MP q, r; mp_divmod(a, b, p, &q, &r); a = b; b = r; }
    return a.d >= 0 ? mp_monic(a, p) : a;
}
static void mp_bezout(MP a, MP b, int64_t p, MP *s, MP *t) {   /* s a + t b = 1 for coprime a, b */
    MP r0 = a, r1 = b, s0 = mp_alloc(0), s1 = mp_alloc(-1), t0 = mp_alloc(-1), t1 = mp_alloc(0);
    s0.c[0] = 1; t1.c[0] = 1;
    while (r1.d >= 0) {
        MP q, r; mp_divmod(r0, r1, p, &q, &r);
        MP s2 = mp_sub(s0, mp_mul(q, s1, p), p), t2 = mp_sub(t0, mp_mul(q, t1, p), p);
        r0 = r1; r1 = r; s0 = s1; s1 = s2; t0 = t1; t1 = t2;
    }
    if (r0.d != 0) nm_fail("internal: Bezout of factors that are not coprime");
    int64_t il = inv_mod(r0.c[0], p);
    *s = mp_alloc(s0.d); for (int i = 0; i <= s0.d; i++) s->c[i] = s0.c[i] * il % p; *s = mp_trim(*s);
    *t = mp_alloc(t0.d); for (int i = 0; i <= t0.d; i++) t->c[i] = t0.c[i] * il % p; *t = mp_trim(*t);
}
static int mp_squarefree(MP f, int64_t p) {
    MP der = mp_alloc(f.d - 1);
    for (int j = 1; j <= f.d; j++) der.c[j - 1] = f.c[j] * j % p;
    der = mp_trim(der);
    return der.d >= 0 && mp_gcd(f, der, p).d == 0;
}

/* Berlekamp: the irreducible factors of a square-free monic g modulo p */
static int berlekamp(MP g, int64_t p, MP *out) {
    int n = g.d;
    int64_t *Qm = arena_alloc((size_t)n * n * sizeof(int64_t));
    MP xp, x = mp_alloc(1), r = mp_alloc(0), q, rem;
    x.c[1] = 1; r.c[0] = 1;
    {
        MP base = x, acc = mp_alloc(0); acc.c[0] = 1;
        for (int64_t e = p; e; ) {
            if (e & 1) { mp_divmod(mp_mul(acc, base, p), g, p, &q, &rem); acc = rem; }
            e >>= 1;
            if (e) { mp_divmod(mp_mul(base, base, p), g, p, &q, &rem); base = rem; }
        }
        xp = acc;
    }
    for (int i = 0; i < n; i++) {                   /* row i: x^(i p) mod g, minus the identity */
        for (int j = 0; j < n; j++) Qm[i * n + j] = j <= r.d ? r.c[j] : 0;
        Qm[i * n + i] = md(Qm[i * n + i] - 1, p);
        mp_divmod(mp_mul(r, xp, p), g, p, &q, &rem); r = rem;
    }
    int64_t *A = arena_alloc((size_t)n * n * sizeof(int64_t));
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[j * n + i] = Qm[i * n + j];
    int *piv = arena_alloc((size_t)(n + 1) * sizeof(int)), rank = 0;
    for (int col = 0; col < n && rank < n; col++) {
        int pr = -1;
        for (int i = rank; i < n; i++) if (A[i * n + col]) { pr = i; break; }
        if (pr < 0) continue;
        for (int j = 0; j < n; j++) { int64_t t = A[pr * n + j]; A[pr * n + j] = A[rank * n + j]; A[rank * n + j] = t; }
        int64_t il = inv_mod(A[rank * n + col], p);
        for (int j = 0; j < n; j++) A[rank * n + j] = A[rank * n + j] * il % p;
        for (int i = 0; i < n; i++) if (i != rank && A[i * n + col]) {
            int64_t f = A[i * n + col];
            for (int j = 0; j < n; j++) A[i * n + j] = md(A[i * n + j] - f * A[rank * n + j], p);
        }
        piv[rank++] = col;
    }
    int k = n - rank;
    MP *basis = arena_alloc((size_t)(k + 1) * sizeof(MP));
    int nb = 0, *isp = arena_alloc((size_t)(n + 1) * sizeof(int));
    for (int j = 0; j < n; j++) isp[j] = -1;
    for (int i = 0; i < rank; i++) isp[piv[i]] = i;
    for (int f = 0; f < n; f++) {
        if (isp[f] >= 0) continue;
        MP v = mp_alloc(n - 1);
        v.c[f] = 1;
        for (int i = 0; i < rank; i++) v.c[piv[i]] = md(-A[i * n + f], p);
        basis[nb++] = mp_trim(v);
    }
    int nf = 1;
    out[0] = g;
    for (int b = 0; b < nb && nf < k; b++) {
        if (basis[b].d <= 0) continue;
        for (int64_t s = 0; s < p && nf < k; s++) {
            MP vs = mp_copy(basis[b]);
            vs.c[0] = md(vs.c[0] - s, p); vs = mp_trim(vs);
            for (int i = 0; i < nf && nf < k; i++) {
                if (out[i].d <= 1) continue;
                MP w = mp_gcd(out[i], vs, p);
                if (w.d > 0 && w.d < out[i].d) {
                    MP qq, rr; mp_divmod(out[i], w, p, &qq, &rr);
                    out[i] = w; out[nf++] = mp_monic(qq, p);
                }
            }
        }
    }
    if (nf != k) nm_fail("internal: Berlekamp's split is incomplete");
    return k;
}

/* ---------------- Hensel's lifting ---------------- */

static ZP zp_from_mp(MP a) { ZP z = zp_alloc(a.d); for (int i = 0; i <= a.d; i++) z.c[i] = zi(a.c[i]); return z; }
static ZP zp_mod(ZP a, Z M) { ZP c = zp_alloc(a.d); for (int i = 0; i <= a.d; i++) c.c[i] = z_mod(a.c[i], M); return zp_trim(c); }
static ZP zp_addscaled(ZP a, MP b, Z P) {          /* a + P b */
    ZP c = zp_alloc(a.d > b.d ? a.d : b.d);
    for (int i = 0; i <= c.d; i++) c.c[i] = z_add(i <= a.d ? a.c[i] : z_zero(), i <= b.d ? z_mul(P, zi(b.c[i])) : z_zero());
    return zp_trim(c);
}
/* T = g h (mod p), g monic, to T = g h (mod p^k) */
static void hensel(ZP T, ZP *g, ZP *h, int64_t p, int k) {
    MP gp = mp_from(*g, p), hp = mp_from(*h, p), s, t;
    mp_bezout(gp, hp, p, &s, &t);
    Z P = zi(p);
    for (int j = 1; j < k; j++) {
        ZP diff = zp_sub(T, zp_mul(*g, *h));
        ZP e = zp_alloc(diff.d);
        for (int i = 0; i <= diff.d; i++) {
            Z q, r; z_divmod(diff.c[i], P, &q, &r);
            if (r.s) nm_fail("internal: Hensel step not divisible");
            e.c[i] = q;
        }
        MP em = mp_from(zp_trim(e), p), q, tau;
        mp_divmod(mp_mul(t, em, p), gp, p, &q, &tau);
        MP sigma = mp_add(mp_mul(s, em, p), mp_mul(q, hp, p), p);
        *g = zp_addscaled(*g, tau, P);
        *h = zp_addscaled(*h, sigma, P);
        P = z_mul(P, zi(p));
    }
    *g = zp_mod(*g, P); *h = zp_mod(*h, P);
}

static Z sym(Z a, Z M) {                           /* the representative in (-M/2, M/2] */
    a = z_mod(a, M);
    if (z_cmp(z_add(a, a), M) > 0) a = z_sub(a, M);
    return a;
}

static Poly zp_to_poly(ZP f, const char *var) {
    Poly p; p.deg = f.d; p.var = var;
    p.c = arena_alloc((size_t)(f.d + 1) * sizeof(Q));
    for (int i = 0; i <= f.d; i++) p.c[i] = q_from_z(f.c[i]);
    return p;
}

/* ---------------- Zassenhaus: the irreducible factors of a primitive square-free f ---------------- */

#define MAXR 40

static int zassenhaus(ZP f, ZP *out, int max) {
    int d = f.d;
    if (d <= 1) { out[0] = f; return 1; }
    Z lc = f.c[d];
    static const int primes[] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
        101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
        211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 0};
    MP *mf = arena_alloc((size_t)(d + 1) * sizeof(MP)), *best = arena_alloc((size_t)(d + 1) * sizeof(MP));
    int r = 0, tried = 0;
    int64_t p = 0;
    for (int i = 0; primes[i] && tried < 5; i++) {    /* of a few good primes, the one with fewest factors */
        int64_t q = primes[i];
        if (zmodp(lc, q) == 0) continue;
        MP fq = mp_from(f, q);
        if (!mp_squarefree(fq, q)) continue;
        tried++;
        int rq = berlekamp(mp_monic(fq, q), q, mf);
        if (!p || rq < r) { r = rq; p = q; MP *t = best; best = mf; mf = t; }
        if (r == 1) break;
    }
    if (!p) nm_fail("internal: no prime below 300 keeps the polynomial square-free");
    mf = best;
    if (nm_show_work) nm_work("modulo %lld: %d factor%s", (long long)p, r, r > 1 ? "s" : "");
    if (r == 1) { out[0] = f; return 1; }
    if (r > MAXR) nm_fail("%d factors modulo %lld: recombination by lattice reduction (van Hoeij) is not built yet", r, (long long)p);
    /* Mignotte: any factor's coefficients are below 2^d |lc| sum |a_i|; lift past twice that */
    Z B = zi(0);
    for (int i = 0; i <= d; i++) B = z_add(B, z_abs(f.c[i]));
    B = z_mul(z_mul(B, z_abs(lc)), z_pow(zi(2), (unsigned)d + 1));
    int k = 1; Z M = zi(p);
    while (z_cmp(M, B) <= 0) { M = z_mul(M, zi(p)); k++; }
    if (nm_show_work) nm_work("lifted to %lld^%d by Hensel's method (past Mignotte's bound)", (long long)p, k);
    ZP *F = arena_alloc((size_t)(r + 1) * sizeof(ZP));
    ZP rest = f;
    for (int i = 0; i < r - 1; i++) {
        MP others = mp_alloc(0); others.c[0] = zmodp(lc, p);
        for (int j = i + 1; j < r; j++) others = mp_mul(others, mf[j], p);
        ZP g = zp_from_mp(mf[i]), h = zp_from_mp(others);
        hensel(rest, &g, &h, p, k);
        F[i] = g; rest = h;
    }
    {
        Z il; if (!z_invmod(lc, M, &il)) nm_fail("internal: leading coefficient not invertible");
        ZP last = zp_alloc(rest.d);
        for (int i = 0; i <= rest.d; i++) last.c[i] = z_mod(z_mul(rest.c[i], il), M);
        F[r - 1] = zp_trim(last);
    }
    /* recombine: subsets by size; a candidate that divides exactly is a true factor */
    int n = 0, used[MAXR] = {0}, left = r;
    ZP g = f;
    long trials = 0;
    for (int s = 1; 2 * s <= left; s++) {
        int idx[MAXR];
        for (int i = 0; i < s; i++) idx[i] = i;
        for (;;) {
            int ok = 1;
            for (int i = 0; i < s; i++) if (used[idx[i]]) ok = 0;
            if (ok) {
                trials++;
                Z glc = g.c[g.d], g0 = g.c[0];
                Z t0 = glc;                           /* constant term first: it must divide lc(g) g(0) */
                for (int i = 0; i < s; i++) t0 = z_mod(z_mul(t0, F[idx[i]].d >= 0 ? F[idx[i]].c[0] : z_zero()), M);
                t0 = sym(t0, M);
                if (g0.s && (t0.s == 0 || z_mod(z_mul(glc, g0), t0).s != 0)) ok = 0;
                if (ok) {
                    ZP cand = zp_alloc(0); cand.c[0] = glc;
                    for (int i = 0; i < s; i++) cand = zp_mod(zp_mul(cand, F[idx[i]]), M);
                    for (int i = 0; i <= cand.d; i++) cand.c[i] = sym(cand.c[i], M);
                    ZP cz = zp_primitive(cand), qz;
                    if (cz.d > 0 && zp_divexact(g, cz, &qz)) {
                        if (n >= max) nm_fail("too many factors");
                        out[n++] = cz;
                        g = qz;
                        for (int i = 0; i < s; i++) used[idx[i]] = 1;
                        left -= s;
                        if (2 * s > left) break;
                    }
                }
            }
            int i = s - 1;                            /* the next subset of size s */
            while (i >= 0 && idx[i] == r - s + i) i--;
            if (i < 0) break;
            idx[i]++;
            for (int j = i + 1; j < s; j++) idx[j] = idx[j - 1] + 1;
        }
    }
    if (g.d > 0) { if (n >= max) nm_fail("too many factors"); out[n++] = zp_primitive(g); }
    if (nm_show_work) nm_work("recombined: %d true factor%s after %ld trial products", n, n > 1 ? "s" : "", trials);
    return n;
}

/* ---------------- the whole factorization ---------------- */

static void add_factor(PolyFactors *R, ZP f, int e, const char *var) {
    if (R->n >= NM_MAXPF) nm_fail("too many factors");
    R->f[R->n] = zp_to_poly(f, var); R->e[R->n++] = e;
}

PolyFactors poly_factor(Poly p) {
    PolyFactors R; memset(&R, 0, sizeof R);
    while (p.deg >= 0 && !q_sign(p.c[p.deg])) p.deg--;
    if (p.deg < 0) nm_fail("the zero polynomial has no factors");
    const char *var = p.var ? p.var : "x";
    Z l = zi(1);                                      /* p = unit * F, F whole and primitive */
    for (int i = 0; i <= p.deg; i++) { Z g = z_gcd(l, p.c[i].den), q, r; z_divmod(z_mul(l, p.c[i].den), g, &q, &r); l = q; }
    ZP F = zp_alloc(p.deg);
    for (int i = 0; i <= p.deg; i++) { Z q, r; z_divmod(z_mul(p.c[i].num, l), p.c[i].den, &q, &r); F.c[i] = q; }
    ZP Fp = zp_primitive(F);
    { Z q, r; z_divmod(F.c[F.d], Fp.c[Fp.d], &q, &r); R.unit = q_make(q, l); }
    F = Fp;
    int low = 0;
    while (!F.c[low].s) low++;
    if (low) {                                        /* x^low */
        ZP x = zp_alloc(1); x.c[1] = zi(1);
        add_factor(&R, x, low, var);
        ZP G = zp_alloc(F.d - low);
        for (int i = 0; i <= G.d; i++) G.c[i] = F.c[i + low];
        F = G;
    }
    if (F.d >= 1) {
        ZP parts[NM_MAXPF];
        ZP Fd = zp_deriv(F), g0 = zp_gcd(F, Fd);
        if (g0.d == 0) {
            int np = zassenhaus(F, parts, NM_MAXPF);
            for (int j = 0; j < np; j++) add_factor(&R, parts[j], 1, var);
        } else {                                      /* Yun: F = prod a_i^i, a_i square-free and coprime */
            ZP b, c, dd;
            if (!zp_divexact(F, g0, &b) || !zp_divexact(Fd, g0, &c)) nm_fail("internal: Yun's division");
            dd = zp_sub(c, zp_deriv(b));
            for (int i = 1; b.d > 0; i++) {
                ZP ai = dd.d < 0 ? zp_primitive(b) : zp_gcd(b, dd);
                if (ai.d > 0) {
                    if (nm_show_work) nm_work("the part to the power %d: %s", i, p_to_str(zp_to_poly(ai, var)));
                    int np = zassenhaus(ai, parts, NM_MAXPF);
                    for (int j = 0; j < np; j++) add_factor(&R, parts[j], i, var);
                }
                ZP nb, nc;
                if (!zp_divexact(b, ai, &nb) || !zp_divexact(dd, ai, &nc)) nm_fail("internal: Yun's division");
                b = nb;
                dd = zp_sub(nc, zp_deriv(b));
            }
        }
    }
    for (int i = 1; i < R.n; i++)                    /* a fixed order: by degree, then as written */
        for (int j = i; j > 0; j--) {
            int a = R.f[j - 1].deg, bd = R.f[j].deg;
            if (a < bd || (a == bd && strcmp(p_to_str(R.f[j - 1]), p_to_str(R.f[j])) <= 0)) break;
            Poly t = R.f[j]; R.f[j] = R.f[j - 1]; R.f[j - 1] = t;
            int te = R.e[j]; R.e[j] = R.e[j - 1]; R.e[j - 1] = te;
        }
    Poly back = p_const(R.unit);                     /* multiplied back */
    for (int i = 0; i < R.n; i++) back = p_mul(back, p_pow(R.f[i], (unsigned)R.e[i]));
    Poly diffp = p_sub(back, p);
    for (int i = 0; i <= diffp.deg; i++) if (q_sign(diffp.c[i])) nm_fail("internal: the factors do not multiply back");
    return R;
}
