/* Series with exact rational coefficients, worked "in speciebus" the way decimals are worked in numbers:
 * c[0] + c[1] x + ... + c[n-1] x^(n-1) + O(x^n). Every operation carries the precision n it can vouch for. */
#include "nm.h"

#include <stdio.h>
#include <string.h>

static Q q0(void) { return q_from_z(z_zero()); }
static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }

static Ser s_alloc(int n) {
    Ser s;
    s.n = n;
    s.c = arena_alloc((size_t)(n > 0 ? n : 1) * sizeof(Q));
    for (int i = 0; i < n; i++) s.c[i] = q0();
    return s;
}

static int mini(int a, int b) { return a < b ? a : b; }

Ser s_const(Q a, int n) { Ser s = s_alloc(n); if (n) s.c[0] = a; return s; }
Ser s_var(int n) { Ser s = s_alloc(n); if (n > 1) s.c[1] = qi(1); return s; }

Ser s_from_poly(Poly p, int n) {
    Ser s = s_alloc(n);
    for (int i = 0; i <= p.deg && i < n; i++) s.c[i] = p.c[i];
    return s;
}

Ser s_trunc(Ser a, int n) { if (n < a.n) a.n = n; return a; }

Ser s_add(Ser a, Ser b) {
    Ser s = s_alloc(mini(a.n, b.n));
    for (int i = 0; i < s.n; i++) s.c[i] = q_add(a.c[i], b.c[i]);
    return s;
}

Ser s_scale(Ser a, Q k) {
    Ser s = s_alloc(a.n);
    for (int i = 0; i < a.n; i++) s.c[i] = q_mul(a.c[i], k);
    return s;
}

Ser s_sub(Ser a, Ser b) { return s_add(a, s_scale(b, qi(-1))); }

Ser s_mul(Ser a, Ser b) {
    Ser s = s_alloc(mini(a.n, b.n));
    for (int i = 0; i < a.n && i < s.n; i++) {
        if (q_sign(a.c[i]) == 0) continue;
        for (int j = 0; i + j < s.n; j++)
            if (q_sign(b.c[j])) s.c[i + j] = q_add(s.c[i + j], q_mul(a.c[i], b.c[j]));
    }
    return s;
}

/* 1/b term by term, as Newton divides a by 1 + x x: each new term cancels what the previous ones left */
static Ser s_inv(Ser b) {
    if (b.n == 0) return b;
    if (q_sign(b.c[0]) == 0) nm_fail("division by a series that vanishes at 0 (negative powers are not in this version)");
    Ser s = s_alloc(b.n);
    Q inv0 = q_div(qi(1), b.c[0]);
    s.c[0] = inv0;
    for (int k = 1; k < b.n; k++) {
        Q acc = q0();
        for (int j = 1; j <= k; j++) if (q_sign(b.c[j])) acc = q_add(acc, q_mul(b.c[j], s.c[k - j]));
        s.c[k] = q_neg(q_mul(acc, inv0));
    }
    return s;
}

Ser s_div(Ser a, Ser b) { return s_mul(a, s_inv(b)); }

Ser s_pow_int(Ser a, int64_t e) {
    if (e < 0) return s_inv(s_pow_int(a, -e));
    Ser r = s_const(qi(1), a.n);
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
    if (q_sign(a.c[0]) == 0) nm_fail("a fractional power of a series that vanishes at 0 (x^(1/2) and the like come later)");
    int64_t p, q;
    if (!z_fits_i64(alpha.num, &p) || !z_fits_i64(alpha.den, &q) || q > 1000) nm_fail("exponent too large");
    Q r0;
    if (!q_root_exact(q_pow(a.c[0], p), q, &r0))
        nm_fail("the series starts with %s, whose power %s is not rational (irrational coefficients come later)",
                q_to_str(a.c[0]), q_to_str(alpha));
    Ser b = s_alloc(a.n);
    b.c[0] = r0;
    Q a0inv = q_div(qi(1), a.c[0]);
    Q ap1 = q_add(alpha, qi(1));
    for (int k = 1; k < a.n; k++) {
        Q acc = q0();
        for (int j = 1; j <= k; j++) {
            if (q_sign(a.c[j]) == 0) continue;
            Q f = q_sub(q_mul(ap1, qi(j)), qi(k));
            acc = q_add(acc, q_mul(f, q_mul(a.c[j], b.c[k - j])));
        }
        b.c[k] = q_mul(acc, q_mul(a0inv, q_div(qi(1), qi(k))));
    }
    /* the direct operation: b^q must give back a^p */
    Ser back = s_pow_int(b, q), want = s_pow_int(a, p);
    for (int i = 0; i < back.n && i < want.n; i++)
        if (q_sign(q_sub(back.c[i], want.c[i]))) nm_fail("internal check failed in a fractional power (vitiose)");
    return b;
}

Ser s_deriv(Ser a) {
    Ser s = s_alloc(a.n > 0 ? a.n - 1 : 0);
    for (int i = 0; i < s.n; i++) s.c[i] = q_mul(a.c[i + 1], qi(i + 1));
    return s;
}

/* the fluent that vanishes at 0: Newton's first rule, term by term */
Ser s_integ(Ser a) {
    Ser s = s_alloc(a.n + 1);
    for (int i = 0; i < a.n; i++) s.c[i + 1] = q_div(a.c[i], qi(i + 1));
    return s;
}

Ser s_compose(Ser f, Ser g) {
    if (g.n && q_sign(g.c[0])) nm_fail("substitution needs a series that vanishes at 0");
    int n = mini(f.n, g.n);
    Ser r = s_const(q0(), n);
    for (int i = f.n - 1; i >= 0; i--) r = s_add(s_mul(r, s_trunc(g, n)), s_const(f.c[i], n));
    return r;
}

Ser s_persist(Ser a) {
    Q *c = perm_alloc((size_t)(a.n > 0 ? a.n : 1) * sizeof(Q));
    for (int i = 0; i < a.n; i++) c[i] = q_persist(a.c[i]);
    a.c = c;
    return a;
}

/* one term as a mathematician writes it: 5x^8/128, x^3/3, 2x, 3/4 */
char *term_str(Q c, const char *var, int k, int first) {
    int neg = q_sign(c) < 0;
    Q a = neg ? q_neg(c) : c;
    char *num = z_to_str(a.num), *den = z_to_str(a.den);
    char *s = arena_alloc(strlen(num) + strlen(den) + strlen(var ? var : "") + 40), *o = s;
    if (first) o += sprintf(o, neg ? "-" : "");
    else o += sprintf(o, neg ? " - " : " + ");
    if (k == 0) { o += sprintf(o, "%s", num); if (!q_is_int(a)) sprintf(o, "/%s", den); return s; }
    if (!z_is_one(a.num)) o += sprintf(o, "%s", num);
    o += sprintf(o, "%s", var);
    if (k > 1) o += sprintf(o, "^%d", k);
    if (!q_is_int(a)) sprintf(o, "/%s", den);
    return s;
}

char *s_to_str(Ser a, const char *var, int last) {
    int upto = a.n - 1 < last ? a.n - 1 : last;
    size_t cap = 64;
    for (int i = 0; i <= upto; i++) cap += strlen(q_to_str(a.c[i])) + strlen(var) + 24;
    char *s = arena_alloc(cap), *o = s;
    *o = 0;
    int first = 1;
    for (int i = 0; i <= upto; i++) {
        if (q_sign(a.c[i]) == 0) continue;
        o += sprintf(o, "%s", term_str(a.c[i], var, i, first));
        first = 0;
    }
    sprintf(o, "%sO(%s^%d)", first ? "" : " + ", var, upto + 1);
    return s;
}
