/* Rationals (reduced fractions) and polynomials in one letter with rational coefficients. */
#include "nm.h"

#include <stdio.h>
#include <string.h>

/* ---------------- rationals ---------------- */

Q q_from_z(Z a) { Q r = {a, z_from_i64(1)}; return r; }

Q q_make(Z num, Z den) {
    if (den.s == 0) nm_fail("division by zero");
    if (den.s < 0) { num = z_neg(num); den = z_neg(den); }
    Z g = z_gcd(num, den), rem;
    if (!z_is_one(g) && num.s) { z_divmod(num, g, &num, &rem); z_divmod(den, g, &den, &rem); }
    if (num.s == 0) den = z_from_i64(1);
    Q r = {num, den};
    return r;
}

Q q_add(Q a, Q b) { return q_make(z_add(z_mul(a.num, b.den), z_mul(b.num, a.den)), z_mul(a.den, b.den)); }
Q q_sub(Q a, Q b) { return q_add(a, q_neg(b)); }
Q q_mul(Q a, Q b) { return q_make(z_mul(a.num, b.num), z_mul(a.den, b.den)); }
Q q_div(Q a, Q b) {
    if (b.num.s == 0) nm_fail("division by zero");
    return q_make(z_mul(a.num, b.den), z_mul(a.den, b.num));
}
Q q_neg(Q a) { a.num = z_neg(a.num); return a; }
int q_sign(Q a) { return a.num.s; }
int q_is_int(Q a) { return z_is_one(a.den); }

Q q_pow(Q a, int64_t e) {
    if (e < 0) {
        if (a.num.s == 0) nm_fail("division by zero");
        Q inv = q_make(a.den, a.num);
        return q_pow(inv, -e);
    }
    if (e > 1000000) nm_fail("exponent too large");
    Q r = {z_pow(a.num, (unsigned)e), z_pow(a.den, (unsigned)e)};
    return r;
}

char *q_to_str(Q a) {
    if (q_is_int(a)) return z_to_str(a.num);
    char *n = z_to_str(a.num), *d = z_to_str(a.den);
    char *s = arena_alloc(strlen(n) + strlen(d) + 2);
    sprintf(s, "%s/%s", n, d);
    return s;
}

Q q_persist(Q a) { a.num = z_persist(a.num); a.den = z_persist(a.den); return a; }

/* ---------------- polynomials ---------------- */

static Poly p_alloc(int deg, const char *var) {
    Poly p;
    p.deg = deg; p.var = var;
    p.c = arena_alloc((size_t)(deg + 1 > 0 ? deg + 1 : 1) * sizeof(Q));
    for (int i = 0; i <= deg; i++) p.c[i] = q_from_z(z_zero());
    return p;
}

static Poly p_trim(Poly p) {
    while (p.deg >= 0 && q_sign(p.c[p.deg]) == 0) p.deg--;
    if (p.deg <= 0) p.var = NULL;
    return p;
}

static const char *join_var(Poly a, Poly b) {
    if (a.var && b.var && strcmp(a.var, b.var) != 0)
        nm_fail("two different unknowns, %s and %s: only one letter per equation in this version", a.var, b.var);
    return a.var ? a.var : b.var;
}

Poly p_const(Q a) { Poly p = p_alloc(0, NULL); p.c[0] = a; return p_trim(p); }

Poly p_var(const char *name) {
    Poly p = p_alloc(1, name);
    p.c[1] = q_from_z(z_from_i64(1));
    return p;
}

Poly p_add(Poly a, Poly b) {
    Poly r = p_alloc(a.deg > b.deg ? a.deg : b.deg, join_var(a, b));
    for (int i = 0; i <= r.deg; i++) {
        Q x = i <= a.deg ? a.c[i] : q_from_z(z_zero());
        Q y = i <= b.deg ? b.c[i] : q_from_z(z_zero());
        r.c[i] = q_add(x, y);
    }
    return p_trim(r);
}

Poly p_scale(Poly a, Q s) {
    Poly r = p_alloc(a.deg, a.var);
    for (int i = 0; i <= a.deg; i++) r.c[i] = q_mul(a.c[i], s);
    return p_trim(r);
}

Poly p_sub(Poly a, Poly b) { return p_add(a, p_scale(b, q_from_z(z_from_i64(-1)))); }

Poly p_mul(Poly a, Poly b) {
    if (a.deg < 0 || b.deg < 0) return p_alloc(-1, NULL);
    Poly r = p_alloc(a.deg + b.deg, join_var(a, b));
    for (int i = 0; i <= a.deg; i++)
        for (int j = 0; j <= b.deg; j++) r.c[i + j] = q_add(r.c[i + j], q_mul(a.c[i], b.c[j]));
    return p_trim(r);
}

Poly p_pow(Poly a, unsigned e) {
    if ((int64_t)a.deg * e > 100000) nm_fail("polynomial degree too large");
    Poly r = p_const(q_from_z(z_from_i64(1)));
    for (unsigned i = 0; i < e; i++) r = p_mul(r, a);
    return r;
}

char *p_to_str(Poly p) {
    if (p.deg < 0) return q_to_str(q_from_z(z_zero()));
    size_t cap = 64;
    for (int i = 0; i <= p.deg; i++) cap += strlen(q_to_str(p.c[i])) + 32;
    char *s = arena_alloc(cap), *o = s;
    *o = 0;
    int first = 1;
    for (int i = p.deg; i >= 0; i--) {
        if (q_sign(p.c[i]) == 0) continue;
        o += sprintf(o, "%s", term_str(p.c[i], p.var, i, first));
        first = 0;
    }
    return s;
}
