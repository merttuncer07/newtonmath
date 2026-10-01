/* Whole numbers: divisors, residues, primes and factors.
 *
 * Newton's notebook counts divisors ("every quantity hath one divisor more than it hath aliquote parts") and
 * copies Descartes' rule for amicable numbers. His standard decides the verdicts here: a prime is "proved" only
 * when a proof was carried out (a certain test below 3.3e24, or Pocklington's certificate from the factors of
 * n - 1, themselves proved); otherwise it is "probable" (the BPSW test, with no known counterexample), and says
 * so. Every factorization is multiplied back. */
#include "nm.h"

#include <string.h>

static Z zi(int64_t v) { return z_from_i64(v); }

Z z_mod(Z a, Z m) {                               /* the residue in [0, |m|) */
    Z q, r;
    z_divmod(a, m, &q, &r);
    if (r.s < 0) r = z_add(r, z_abs(m));
    return r;
}

Z z_powmod(Z b, Z e, Z m) {
    if (e.s < 0) nm_fail("a negative power modulo m needs an inverse");
    Z r = z_mod(zi(1), m);
    b = z_mod(b, m);
    /* bits of e, from the top */
    Z bits[4096]; int nb = 0;
    Z two = zi(2), q, rem;
    while (e.s) { z_divmod(e, two, &q, &rem); bits[nb++] = rem; e = q; if (nb == 4096) nm_fail("exponent too large"); }
    for (int i = nb - 1; i >= 0; i--) {
        r = z_mod(z_mul(r, r), m);
        if (bits[i].s) r = z_mod(z_mul(r, b), m);
    }
    return r;
}

Z z_xgcd(Z a, Z b, Z *s, Z *t) {                  /* g = s a + t b */
    Z s0 = zi(1), s1 = zi(0), t0 = zi(0), t1 = zi(1);
    while (b.s) {
        Z q, r;
        z_divmod(a, b, &q, &r);
        Z s2 = z_sub(s0, z_mul(q, s1)), t2 = z_sub(t0, z_mul(q, t1));
        a = b; b = r; s0 = s1; s1 = s2; t0 = t1; t1 = t2;
    }
    if (a.s < 0) { a = z_neg(a); s0 = z_neg(s0); t0 = z_neg(t0); }
    *s = s0; *t = t0;
    return a;
}

int z_invmod(Z a, Z m, Z *out) {
    Z s, t, g = z_xgcd(z_mod(a, m), m, &s, &t);
    if (!z_is_one(g)) return 0;
    *out = z_mod(s, m);
    return 1;
}

/* x = r_i mod m_i for all i; returns 0 if they contradict each other */
int z_crt(const Z *r, const Z *m, int k, Z *x, Z *M) {
    Z X = zi(0), MM = zi(1);
    for (int i = 0; i < k; i++) {
        Z s, t, g = z_xgcd(MM, m[i], &s, &t), q, rem;
        z_divmod(z_sub(r[i], X), g, &q, &rem);
        if (rem.s) return 0;
        Z l, dummy;
        z_divmod(m[i], g, &l, &dummy);
        Z step = z_mod(z_mul(q, s), l);
        X = z_add(X, z_mul(MM, step));
        MM = z_mul(MM, l);
        X = z_mod(X, MM);
    }
    *x = X; *M = MM;
    return 1;
}

/* ---------------- primes ---------------- */

static const int SMALL[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41};

static int strong_probable(Z n, Z a) {           /* Miller-Rabin to base a */
    Z one = zi(1), nm1 = z_sub(n, one), d = nm1, q, r;
    int s = 0;
    for (;;) { z_divmod(d, zi(2), &q, &r); if (r.s) break; d = q; s++; }
    Z x = z_powmod(a, d, n);
    if (z_is_one(x) || !z_cmp(x, nm1)) return 1;
    for (int i = 1; i < s; i++) {
        x = z_mod(z_mul(x, x), n);
        if (!z_cmp(x, nm1)) return 1;
        if (z_is_one(x)) return 0;
    }
    return 0;
}

static int jacobi(Z a, Z n) {
    a = z_mod(a, n);
    int t = 1;
    Z q, r;
    while (a.s) {
        for (;;) {
            z_divmod(a, zi(2), &q, &r);
            if (r.s) break;
            a = q;
            Z m8 = z_mod(n, zi(8));
            int64_t v; z_fits_i64(m8, &v);
            if (v == 3 || v == 5) t = -t;
        }
        Z tmp = a; a = n; n = tmp;
        int64_t va, vn; z_fits_i64(z_mod(a, zi(4)), &va); z_fits_i64(z_mod(n, zi(4)), &vn);
        if (va == 3 && vn == 3) t = -t;
        a = z_mod(a, n);
    }
    return z_is_one(n) ? t : 0;
}

static int strong_lucas(Z n) {                    /* Selfridge's parameters, the strong Lucas test */
    int64_t D = 5;
    for (int i = 0; i < 1000; i++) {
        int j = jacobi(zi(D), n);
        if (j == -1) break;
        if (j == 0 && z_cmp(z_abs(zi(D)), n) != 0) return 0;
        D = D > 0 ? -(D + 2) : -(D - 2);
        if (i == 20) {                            /* a perfect square has no such D */
            Z s = z_iroot(n, 2);
            if (!z_cmp(z_mul(s, s), n)) return 0;
        }
    }
    Z P = zi(1), Qp = z_from_i64((1 - D) / 4);
    Z d = z_add(n, zi(1)), q, r;
    int s = 0;
    for (;;) { z_divmod(d, zi(2), &q, &r); if (r.s) break; d = q; s++; }
    /* U_d, V_d by doubling along the bits of d */
    Z bits[4096]; int nb = 0;
    Z e = d;
    while (e.s) { z_divmod(e, zi(2), &q, &r); bits[nb++] = r; e = q; }
    Z U = zi(1), V = P, Qk = z_mod(Qp, n);
    Z inv2; z_invmod(zi(2), n, &inv2);
    for (int i = nb - 2; i >= 0; i--) {
        U = z_mod(z_mul(U, V), n);
        V = z_mod(z_sub(z_mul(V, V), z_mul(zi(2), Qk)), n);
        Qk = z_mod(z_mul(Qk, Qk), n);
        if (bits[i].s) {
            Z U2 = z_mod(z_mul(z_add(z_mul(P, U), V), inv2), n);
            Z V2 = z_mod(z_mul(z_add(z_mul(zi(D), U), z_mul(P, V)), inv2), n);
            U = U2; V = V2;
            Qk = z_mod(z_mul(Qk, Qp), n);
        }
    }
    if (!U.s || !V.s) return 1;
    for (int i = 1; i < s; i++) {
        V = z_mod(z_sub(z_mul(V, V), z_mul(zi(2), Qk)), n);
        if (!V.s) return 1;
        Qk = z_mod(z_mul(Qk, Qk), n);
    }
    return 0;
}

static int prove_pocklington(Z n, int depth);

/* 0: composite, 1: proved prime, 2: probable prime */
int z_isprime(Z n) {
    if (n.s <= 0) return 0;
    int64_t v;
    if (z_fits_i64(n, &v) && v < 2) return 0;
    for (int i = 0; i < 13; i++) {
        if (z_fits_i64(n, &v) && v == SMALL[i]) return 1;
        Z q, r; z_divmod(n, zi(SMALL[i]), &q, &r);
        if (!r.s) return 0;
    }
    for (int i = 0; i < 13; i++) if (!strong_probable(n, zi(SMALL[i]))) return 0;
    /* below 3.3 * 10^24 the thirteen bases decide (Sorenson and Webster) */
    if (z_cmp(n, z_mul(z_from_i64(3317044064679887LL), z_from_i64(1000000000LL))) < 0) return 1;
    if (!strong_lucas(n)) return 0;
    return prove_pocklington(n, 0) ? 1 : 2;
}

/* ---------------- factors ---------------- */

static Z rho(Z n) {                               /* a nontrivial divisor of the composite n (Pollard, with Brent's cycle) */
    for (int64_t c = 1; c < 200; c++) {
        Z y = zi(2), x = y, g = zi(1), q = zi(1), ys = y, C = zi(c);
        int64_t r = 1, m = 64;
        do {
            x = y;
            for (int64_t i = 0; i < r; i++) y = z_mod(z_add(z_mul(y, y), C), n);
            int64_t k = 0;
            do {
                ys = y;
                for (int64_t i = 0; i < m && i < r - k; i++) {
                    y = z_mod(z_add(z_mul(y, y), C), n);
                    q = z_mod(z_mul(q, z_abs(z_sub(x, y))), n);
                }
                g = z_gcd(q, n);
                k += m;
            } while (k < r && z_is_one(g));
            r *= 2;
            if (r > (1 << 26)) break;
        } while (z_is_one(g));
        if (!z_cmp(g, n)) {
            do { ys = z_mod(z_add(z_mul(ys, ys), C), n); g = z_gcd(z_abs(z_sub(x, ys)), n); } while (z_is_one(g));
        }
        if (z_cmp(g, n) && !z_is_one(g)) return g;
    }
    nm_fail("the factors of %s could not be separated in reasonable time", z_to_str(n));
    return n;
}

static void add_factor(Factors *F, Z p, int e) {
    for (int i = 0; i < F->n; i++) if (!z_cmp(F->p[i], p)) { F->e[i] += e; return; }
    if (F->n == NM_MAXFAC) nm_fail("too many prime factors");
    int k = F->n++;
    F->p[k] = p; F->e[k] = e;
    F->status[k] = z_isprime(p);
    /* keep them in order */
    while (k > 0 && z_cmp(F->p[k - 1], F->p[k]) > 0) {
        Z tp = F->p[k]; F->p[k] = F->p[k - 1]; F->p[k - 1] = tp;
        int te = F->e[k]; F->e[k] = F->e[k - 1]; F->e[k - 1] = te;
        int ts = F->status[k]; F->status[k] = F->status[k - 1]; F->status[k - 1] = ts;
        k--;
    }
}

static void split(Factors *F, Z n) {
    if (z_is_one(n)) return;
    if (z_isprime(n)) { add_factor(F, n, 1); return; }
    Z s = z_iroot(n, 2);                              /* a square first: rho is slow on them */
    if (!z_cmp(z_mul(s, s), n)) { split(F, s); split(F, s); return; }
    Z d = rho(n), q, r;
    z_divmod(n, d, &q, &r);
    split(F, d); split(F, q);
}

Factors z_factor(Z n) {
    Factors F;
    memset(&F, 0, sizeof F);
    if (n.s == 0) nm_fail("0 has no factorization");
    F.sign = n.s;
    n = z_abs(n);
    Z whole = n;
    for (int64_t p = 2; p < 10000; p += (p == 2 ? 1 : 2)) {   /* small primes by trial */
        int64_t v;
        if (z_fits_i64(n, &v) && p * p > v) break;
        int e = 0;
        for (;;) { Z q, r; z_divmod(n, zi(p), &q, &r); if (r.s) break; n = q; e++; }
        if (e) add_factor(&F, zi(p), e);
    }
    split(&F, n);
    /* the direct operation: the product must give n back */
    Z prod = zi(1);
    for (int i = 0; i < F.n; i++) prod = z_mul(prod, z_pow(F.p[i], (unsigned)F.e[i]));
    if (z_cmp(prod, whole)) nm_fail("internal check failed: the factors do not multiply back (vitiose)");
    return F;
}

/* Pocklington: n is prime if a^(n-1) = 1 and gcd(a^((n-1)/q) - 1, n) = 1 for every prime q of a proved part
 * F of n - 1 with F > sqrt(n). The primes of F are proved the same way. */
static int prove_pocklington(Z n, int depth) {
    if (depth > 20) return 0;
    Z nm1 = z_sub(n, zi(1));
    Z rest = nm1, Fpart = zi(1);
    Z qs[NM_MAXFAC]; int nq = 0;
    for (int64_t p = 2; p < 100000 && nq < NM_MAXFAC; p += (p == 2 ? 1 : 2)) {
        Z q, r; z_divmod(rest, zi(p), &q, &r);
        if (r.s) continue;
        qs[nq++] = zi(p);
        while (!r.s) { rest = q; Fpart = z_mul(Fpart, zi(p)); z_divmod(rest, zi(p), &q, &r); }
    }
    if (z_cmp(z_mul(Fpart, Fpart), n) <= 0 && !z_is_one(rest)) {
        /* the remaining cofactor: prime (proved), or split further */
        if (z_isprime(rest) == 1) { qs[nq++] = rest; Fpart = z_mul(Fpart, rest); rest = zi(1); }
        else {
            Factors G = z_factor(rest);
            for (int i = 0; i < G.n && nq < NM_MAXFAC; i++) {
                if (G.status[i] != 1) continue;
                qs[nq++] = G.p[i];
                Fpart = z_mul(Fpart, z_pow(G.p[i], (unsigned)G.e[i]));
            }
        }
    }
    if (z_cmp(z_mul(Fpart, Fpart), n) <= 0) return 0;
    for (int64_t a = 2; a < 200; a++) {
        if (!z_is_one(z_powmod(zi(a), nm1, n))) return 0;   /* not even a Fermat witness: composite */
        int all = 1;
        for (int i = 0; i < nq && all; i++) {
            Z e, rm; z_divmod(nm1, qs[i], &e, &rm);
            Z g = z_gcd(z_sub(z_powmod(zi(a), e, n), zi(1)), n);
            if (!z_is_one(g)) all = 0;
        }
        if (all) return 1;
    }
    return 0;
}

Z z_sigma(Z n) {                                   /* the sum of all divisors */
    Factors F = z_factor(n);
    Z s = zi(1);
    for (int i = 0; i < F.n; i++) {
        Z num = z_sub(z_pow(F.p[i], (unsigned)(F.e[i] + 1)), zi(1)), q, r;
        z_divmod(num, z_sub(F.p[i], zi(1)), &q, &r);
        s = z_mul(s, q);
    }
    return s;
}

Z z_phi(Z n) {
    Factors F = z_factor(n);
    Z s = zi(1);
    for (int i = 0; i < F.n; i++) s = z_mul(s, z_mul(z_pow(F.p[i], (unsigned)(F.e[i] - 1)), z_sub(F.p[i], zi(1))));
    return s;
}

int z_divisors(Z n, Z *out, int max) {             /* in increasing order */
    Factors F = z_factor(n);
    int k = 1;
    out[0] = zi(1);
    for (int i = 0; i < F.n; i++) {
        int cur = k;
        Z pw = zi(1);
        for (int e = 1; e <= F.e[i]; e++) {
            pw = z_mul(pw, F.p[i]);
            for (int j = 0; j < cur; j++) { if (k == max) nm_fail("too many divisors to list"); out[k++] = z_mul(out[j], pw); }
        }
    }
    for (int i = 1; i < k; i++) for (int j = i; j > 0 && z_cmp(out[j - 1], out[j]) > 0; j--) { Z t = out[j]; out[j] = out[j - 1]; out[j - 1] = t; }
    return k;
}

Z z_nextprime(Z n) {
    Z p = z_add(n, zi(1));
    while (!z_isprime(p)) p = z_add(p, zi(1));
    return p;
}
