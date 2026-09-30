/* Series whose coefficients are quantities in letters, worked "in speciebus" the way decimals are worked in
 * numbers: c[0] + c[1] x + ... + c[n-1] x^(n-1) + O(x^n). Every operation carries the precision n it can vouch
 * for. Division needs the first coefficient to be a single term, as in Newton's examples (a a/(b + x)). */
#include "nm.h"

#include <stdio.h>
#include <string.h>

static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }

static Ser s_alloc(int n) {
    Ser s;
    s.n = n;
    s.c = arena_alloc((size_t)(n > 0 ? n : 1) * sizeof(C));
    for (int i = 0; i < n; i++) s.c[i] = c_zero();
    return s;
}

static int mini(int a, int b) { return a < b ? a : b; }

Ser s_const(C a, int n) { Ser s = s_alloc(n); if (n) s.c[0] = a; return s; }
Ser s_var(int n) { Ser s = s_alloc(n); if (n > 1) s.c[1] = c_const(qi(1)); return s; }

Ser s_trunc(Ser a, int n) { if (n < a.n) a.n = n; return a; }

Ser s_add(Ser a, Ser b) {
    Ser s = s_alloc(mini(a.n, b.n));
    for (int i = 0; i < s.n; i++) s.c[i] = c_add(a.c[i], b.c[i]);
    return s;
}

Ser s_scale(Ser a, C k) {
    Ser s = s_alloc(a.n);
    for (int i = 0; i < a.n; i++) s.c[i] = c_mul(a.c[i], k);
    return s;
}

Ser s_sub(Ser a, Ser b) { return s_add(a, s_scale(b, c_const(qi(-1)))); }

Ser s_mul(Ser a, Ser b) {
    Ser s = s_alloc(mini(a.n, b.n));
    for (int i = 0; i < a.n && i < s.n; i++) {
        if (c_is_zero(a.c[i])) continue;
        for (int j = 0; i + j < s.n; j++)
            if (!c_is_zero(b.c[j])) s.c[i + j] = c_add(s.c[i + j], c_mul(a.c[i], b.c[j]));
    }
    return s;
}

/* 1/b term by term, as Newton divides aa by b + x: each new term cancels what the previous ones left */
static Ser s_inv(Ser b) {
    if (b.n == 0) return b;
    if (c_is_zero(b.c[0])) nm_fail("division by a series that vanishes at 0 (negative powers come with the parallelogram)");
    Ser s = s_alloc(b.n);
    C one = c_const(qi(1));
    s.c[0] = c_div(one, b.c[0]);
    for (int k = 1; k < b.n; k++) {
        C acc = c_zero();
        for (int j = 1; j <= k; j++) if (!c_is_zero(b.c[j])) acc = c_add(acc, c_mul(b.c[j], s.c[k - j]));
        s.c[k] = c_neg(c_mul(acc, s.c[0]));
    }
    return s;
}

Ser s_div(Ser a, Ser b) { return s_mul(a, s_inv(b)); }

Ser s_pow_int(Ser a, int64_t e) {
    if (e < 0) return s_inv(s_pow_int(a, -e));
    Ser r = s_const(c_const(qi(1)), a.n);
    while (e) {
        if (e & 1) r = s_mul(r, a);
        e >>= 1;
        if (e) a = s_mul(a, a);
    }
    return r;
}

int q_root_exact(Q a, int64_t n, Q *out) {
    int neg = q_sign(a) < 0;
    if (neg && n % 2 == 0) return 0;
    Z A = z_abs(a.num), B = a.den;
    Z ra = z_iroot(A, (unsigned)n), rb = z_iroot(B, (unsigned)n);
    if (z_cmp(z_pow(ra, (unsigned)n), A) || z_cmp(z_pow(rb, (unsigned)n), B)) return 0;
    Q r = q_make(ra, rb);
    *out = neg ? q_neg(r) : r;
    return 1;
}

/* a^alpha by the rule that generalises Newton's binomial theorem: each term from the ones before,
 *   b_0 = a_0^alpha,  b_k = (1 / (k a_0)) * sum_{j=1..k} ((alpha + 1) j - k) a_j b_{k-j}.
 * The result is then raised back to the denominator of alpha and compared, as Newton squared his series. */
Ser s_pow_q(Ser a, Q alpha) {
    if (q_is_int(alpha)) {
        int64_t e;
        if (!z_fits_i64(alpha.num, &e) || e > 100000 || e < -100000) nm_fail("exponent too large");
        return s_pow_int(a, e);
    }
    if (a.n == 0) return a;
    if (c_is_zero(a.c[0])) nm_fail("a fractional power of a series that vanishes at 0: write it as a root of an equation (the parallelogram)");
    int64_t p, q;
    if (!z_fits_i64(alpha.num, &p) || !z_fits_i64(alpha.den, &q) || q > 1000) nm_fail("exponent too large");
    C r0;
    if (!c_pow_q(a.c[0], alpha, &r0))
        nm_fail("the series starts with %s, whose power %s is not exact here (irrational numbers come later)",
                c_to_str(a.c[0]), q_to_str(alpha));
    Ser b = s_alloc(a.n);
    b.c[0] = r0;
    C a0inv = c_div(c_const(qi(1)), a.c[0]);
    Q ap1 = q_add(alpha, qi(1));
    for (int k = 1; k < a.n; k++) {
        C acc = c_zero();
        for (int j = 1; j <= k; j++) {
            if (c_is_zero(a.c[j])) continue;
            Q f = q_sub(q_mul(ap1, qi(j)), qi(k));
            acc = c_add(acc, c_scale(c_mul(a.c[j], b.c[k - j]), f));
        }
        b.c[k] = c_scale(c_mul(acc, a0inv), q_div(qi(1), qi(k)));
    }
    /* the direct operation: b^q must give back a^p */
    Ser back = s_pow_int(b, q), want = s_pow_int(a, p);
    for (int i = 0; i < back.n && i < want.n; i++)
        if (!c_equal(back.c[i], want.c[i])) nm_fail("internal check failed in a fractional power (vitiose)");
    return b;
}

Ser s_deriv(Ser a) {
    Ser s = s_alloc(a.n > 0 ? a.n - 1 : 0);
    for (int i = 0; i < s.n; i++) s.c[i] = c_scale(a.c[i + 1], qi(i + 1));
    return s;
}

/* the fluent that vanishes at 0: Newton's first rule, term by term */
Ser s_integ(Ser a) {
    Ser s = s_alloc(a.n + 1);
    for (int i = 0; i < a.n; i++) s.c[i + 1] = c_scale(a.c[i], q_div(qi(1), qi(i + 1)));
    return s;
}

Ser s_compose(Ser f, Ser g) {
    if (g.n && !c_is_zero(g.c[0])) nm_fail("substitution needs a series that vanishes at 0");
    int n = mini(f.n, g.n);
    Ser r = s_const(c_zero(), n);
    for (int i = f.n - 1; i >= 0; i--) r = s_add(s_mul(r, s_trunc(g, n)), s_const(f.c[i], n));
    return r;
}

Ser s_persist(Ser a) {
    C *c = perm_alloc((size_t)(a.n > 0 ? a.n : 1) * sizeof(C));
    for (int i = 0; i < a.n; i++) c[i] = c_persist(a.c[i]);
    a.c = c;
    return a;
}

char *s_to_str(Ser a, const char *var, int last) {
    int upto = a.n - 1 < last ? a.n - 1 : last;
    size_t cap = 64;
    char **parts = arena_alloc((size_t)(upto + 2) * sizeof(char *));
    int first = 1;
    for (int i = 0; i <= upto; i++) {
        parts[i] = NULL;
        if (c_is_zero(a.c[i])) continue;
        parts[i] = c_term_str(a.c[i], var, qi(i), first);
        cap += strlen(parts[i]);
        first = 0;
    }
    char *s = arena_alloc(cap + strlen(var) + 32), *o = s;
    *o = 0;
    for (int i = 0; i <= upto; i++) if (parts[i]) o += sprintf(o, "%s", parts[i]);
    sprintf(o, "%sO(%s^%d)", first ? "" : " + ", var, upto + 1);
    return s;
}
