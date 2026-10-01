/* Finite quotients first; continual division reduces them to lowest terms.
 * The recursive primitive pseudo-remainder sequence uses elimination's UP, with
 * one fewer ordinary letter in each content calculation. No sampling enters gcd.
 * Every gcd is divided into both inputs; every quotient is multiplied back. */
#include "nm.h"

#include <stdio.h>
#include <string.h>

static Q qi(int64_t n) { return q_from_z(z_from_i64(n)); }
static C one(void) { return c_const(qi(1)); }
static int unit(C a) { Q q; return c_const_value(a, &q) && q_cmp_one(q) == 0; }

static void powers(C a, int polynomial) {
    for (int t = 0; t < a.nt; t++) for (int l = 0; l < letter_count(); l++) {
        Q e = a.t[t].e[l];
        if (!q_is_int(e)) nm_fail("fractional exponents in a quotient are not supported");
        if (polynomial && q_sign(e) < 0) nm_fail("polynomial gcd needs nonnegative whole exponents");
    }
}
static int special(C a) {
    for (int l = 0; l < letter_count(); l++)
        if ((letter_is_surd(l) || letter_is_named(l)) && c_uses(a, l)) return 1;
    return 0;
}
static void polynomial(C a) {
    powers(a, 1);
    if (special(a)) nm_fail("polynomial gcd is over rational coefficients and ordinary letters only");
}
static C monic(C a) {
    int k = c_leading_index(a);
    return k < 0 ? a : c_scale(a, q_div(qi(1), a.t[k].k));
}

C poly_exact_div(C a, C b) {
    if (c_is_zero(b)) nm_fail("division by zero polynomial");
    C original = a, q = c_zero();
    Q k;
    if (c_const_value(b, &k)) q = c_scale(a, q_div(qi(1), k));
    else {
        CT lead = b.t[c_leading_index(b)];
        while (!c_is_zero(a)) {
            CT term = a.t[c_leading_index(a)];
            term.k = q_div(term.k, lead.k);
            for (int l = 0; l < NM_MAXL; l++) {
                term.e[l] = q_sub(term.e[l], lead.e[l]);
                if (q_sign(term.e[l]) < 0) nm_fail("internal check failed: polynomial division has a remainder (vitiose)");
            }
            C step = {1, &term};
            q = c_add(q, step);
            a = c_sub(a, c_mul(step, b));
        }
    }
    if (!c_equal(c_mul(q, b), original)) nm_fail("internal check failed: exact polynomial division (vitiose)");
    return q;
}

static C content(UP a) {
    C g = c_zero();
    for (int i = 0; i <= a.deg; i++) {
        g = poly_gcd(g, a.c[i]);
        if (unit(g)) break;
    }
    return g;
}
static UP primitive(UP a, C g) {
    if (a.deg < 0 || unit(g)) return a;
    UP p = up_alloc(a.deg);
    for (int i = 0; i <= a.deg; i++) p.c[i] = poly_exact_div(a.c[i], g);
    return p;
}
static UP pseudo_rem(UP a, UP b) {
    while (a.deg >= b.deg && a.deg >= 0) {
        int sh = a.deg - b.deg;
        UP scaled = up_alloc(a.deg), multiple = up_alloc(a.deg);
        for (int i = 0; i <= a.deg; i++) scaled.c[i] = c_mul(a.c[i], b.c[b.deg]);
        for (int i = 0; i <= b.deg; i++) multiple.c[i + sh] = c_mul(b.c[i], a.c[a.deg]);
        UP r = up_sub(scaled, multiple);
        if (r.deg >= a.deg) nm_fail("internal check failed: pseudo-division did not lower the degree (vitiose)");
        a = r;
    }
    return a;
}
static C up_join(UP a, int v) {
    C p = c_zero(), x = c_letter(v);
    for (int i = a.deg; i >= 0; i--) p = c_add(c_mul(p, x), a.c[i]);
    return p;
}

C poly_gcd(C a, C b) {
    polynomial(a); polynomial(b);
    C g;
    Q q;
    if (c_is_zero(a)) g = monic(b);
    else if (c_is_zero(b)) g = monic(a);
    else if (c_const_value(a, &q) || c_const_value(b, &q)) g = one();
    else if (c_equal(a, b)) g = monic(a);
    else {
        int v = -1;
        for (int l = 0; l < letter_count(); l++)
            if ((c_uses(a, l) || c_uses(b, l)) && (v < 0 || strcmp(letter_name(l), letter_name(v)) > 0)) v = l;
        UP A = up_from(a, v), B = up_from(b, v);
        C ca = content(A), cb = content(B);
        C cg = poly_gcd(ca, cb);
        A = primitive(A, ca); B = primitive(B, cb);
        while (B.deg >= 0) {
            UP r = pseudo_rem(A, B);
            if (r.deg >= 0) r = primitive(r, content(r));
            A = B; B = r;
        }
        g = monic(c_mul(cg, up_join(A, v)));
    }
    if (!c_is_zero(g)) { (void)poly_exact_div(a, g); (void)poly_exact_div(b, g); }
    return g;
}

/* Clear Laurent powers in both sides by the same monomial before polynomial gcd. */
static void clear_negative(C *num, C *den) {
    C shift = one();
    for (int l = 0; l < letter_count(); l++) {
        Q low = qi(0);
        for (int t = 0; t < num->nt; t++) if (q_cmp(num->t[t].e[l], low) < 0) low = num->t[t].e[l];
        for (int t = 0; t < den->nt; t++) if (q_cmp(den->t[t].e[l], low) < 0) low = den->t[t].e[l];
        shift.t[0].e[l] = q_neg(low);
    }
    if (!unit(shift)) { *num = c_mul(*num, shift); *den = c_mul(*den, shift); }
}
static R checked(R r, C num, C den) {
    if (!c_equal(c_mul(r.num, den), c_mul(num, r.den)))
        nm_fail("internal check failed: quotient cross-multiplication (vitiose)");
    return r;
}

R r_from_c(C a) {
    R r = {a, one()};
    if (!special(a)) for (int t = 0; t < a.nt; t++) for (int l = 0; l < letter_count(); l++)
        if (q_is_int(a.t[t].e[l]) && q_sign(a.t[t].e[l]) < 0) return r_make(a, r.den);
    return r;
}
R r_make(C num, C den) {
    if (c_is_zero(den)) nm_fail("division by zero");
    powers(num, 0); powers(den, 0);
    C original_num = num, original_den = den;
    R r;
    if (special(den) && c_has_plain(den)) nm_fail("surds or named numbers mixed with letters in a denominator are not supported");
    if (!special(num) && !special(den)) clear_negative(&num, &den);
    if (!c_has_plain(den)) {
        r.num = c_div(num, den); r.den = one();
        return checked(r, original_num, original_den);
    }
    if (special(num)) nm_fail("a quotient with a letter denominator requires rational coefficients (surds are not supported here)");
    clear_negative(&num, &den);
    C g = poly_gcd(num, den);
    r.num = poly_exact_div(num, g); r.den = poly_exact_div(den, g);
    Q scale = q_div(qi(1), r.den.t[c_leading_index(r.den)].k);
    r.num = c_scale(r.num, scale); r.den = c_scale(r.den, scale);
    return checked(r, original_num, original_den);
}

R r_add(R a, R b) {
    if (unit(a.den) && unit(b.den)) return r_from_c(c_add(a.num, b.num));
    C num = c_add(c_mul(a.num, b.den), c_mul(b.num, a.den));
    return r_make(num, c_mul(a.den, b.den));
}
R r_sub(R a, R b) {
    if (unit(a.den) && unit(b.den)) return r_from_c(c_sub(a.num, b.num));
    C num = c_sub(c_mul(a.num, b.den), c_mul(b.num, a.den));
    return r_make(num, c_mul(a.den, b.den));
}
R r_mul(R a, R b) {
    if (unit(a.den) && unit(b.den)) return r_from_c(c_mul(a.num, b.num));
    return r_make(c_mul(a.num, b.num), c_mul(a.den, b.den));
}
R r_div(R a, R b) {
    if (c_is_zero(b.num)) nm_fail("division by zero");
    return r_make(c_mul(a.num, b.den), c_mul(a.den, b.num));
}
R r_pow_int(R a, int64_t e) {
    uint64_t n = e < 0 ? (uint64_t)(-(e + 1)) + 1 : (uint64_t)e;
    if (e < 0) a = r_make(a.den, a.num);
    R r = r_from_c(one());
    while (n) { if (n & 1) r = r_mul(r, a); n >>= 1; if (n) a = r_mul(a, a); }
    return r;
}
int r_equal(R a, R b) { return c_equal(c_mul(a.num, b.den), c_mul(b.num, a.den)); }
int r_is_zero(R a) { return c_is_zero(a.num); }
R r_persist(R a) { R r = {c_persist(a.num), c_persist(a.den)}; return r; }

char *r_to_str(R a) {
    if (unit(a.den)) return c_to_str(a.num);
    char *n = c_to_str(a.num), *d = c_to_str(a.den);
    /* A coefficient fraction also needs grouping inside an outer quotient. */
    int np = a.num.nt > 1 || strchr(n, '/') != NULL;
    int dp = a.den.nt > 1 || strchr(d, '/') != NULL;
    char *s = arena_alloc(strlen(n) + strlen(d) + 8);
    sprintf(s, "%s%s%s/%s%s%s", np ? "(" : "", n, np ? ")" : "", dp ? "(" : "", d, dp ? ")" : "");
    return s;
}
