/* Quantities in letters: sums of terms  k · a^e1 · b^e2 ...  with rational k and rational exponents, so that
 * a^(1/2) and a^(-1) are written as Newton wrote them, "ob analogiam rei". Division is by a single term only,
 * as in all of Newton's examples (xx/64a, 131x^3/512aa); anything else is refused with the reason. */
#include "nm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------- letters ---------------- */

static char *letter_names[NM_MAXL];
static int nletters;

int letter_index(const char *name, size_t len) {
    for (int i = 0; i < nletters; i++)
        if (strlen(letter_names[i]) == len && !strncmp(letter_names[i], name, len)) return i;
    if (nletters == NM_MAXL) nm_fail("too many different letters (at most %d)", NM_MAXL);
    char *s = perm_alloc(len + 1);
    memcpy(s, name, len); s[len] = 0;
    letter_names[nletters] = s;
    return nletters++;
}

const char *letter_name(int i) { return letter_names[i]; }
int letter_count(void) { return nletters; }

/* ---------------- terms ---------------- */

static Q q0(void) { return q_from_z(z_zero()); }
static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }

static int exps_cmp(const CT *a, const CT *b) {
    for (int i = 0; i < NM_MAXL; i++) {
        int c = q_cmp(a->e[i], b->e[i]);
        if (c) return c;
    }
    return 0;
}

static int cmp_terms(const void *x, const void *y) { return exps_cmp(x, y); }

static CT term_one(void) {
    CT t;
    t.k = qi(1);
    for (int i = 0; i < NM_MAXL; i++) t.e[i] = q0();
    return t;
}

static C c_alloc(int nt) {
    C c;
    c.nt = nt;
    c.t = arena_alloc((size_t)(nt > 0 ? nt : 1) * sizeof(CT));
    return c;
}

static C c_canon(C c) {                       /* sort, merge equal exponents, drop zeros */
    if (c.nt > 1) qsort(c.t, (size_t)c.nt, sizeof(CT), cmp_terms);
    int w = 0;
    for (int i = 0; i < c.nt; i++) {
        if (w > 0 && exps_cmp(&c.t[w - 1], &c.t[i]) == 0) c.t[w - 1].k = q_add(c.t[w - 1].k, c.t[i].k);
        else c.t[w++] = c.t[i];
    }
    int z = 0;
    for (int i = 0; i < w; i++) if (q_sign(c.t[i].k)) c.t[z++] = c.t[i];
    c.nt = z;
    return c;
}

C c_zero(void) { return c_alloc(0); }

C c_const(Q k) {
    if (q_sign(k) == 0) return c_zero();
    C c = c_alloc(1);
    c.t[0] = term_one();
    c.t[0].k = k;
    return c;
}

C c_letter(int idx) {
    C c = c_alloc(1);
    c.t[0] = term_one();
    c.t[0].e[idx] = qi(1);
    return c;
}

int c_is_zero(C a) { return a.nt == 0; }

int c_const_value(C a, Q *out) {
    if (a.nt == 0) { *out = q0(); return 1; }
    if (a.nt > 1) return 0;
    for (int i = 0; i < NM_MAXL; i++) if (q_sign(a.t[0].e[i])) return 0;
    *out = a.t[0].k;
    return 1;
}

int c_uses(C a, int idx) {
    for (int i = 0; i < a.nt; i++) if (q_sign(a.t[i].e[idx])) return 1;
    return 0;
}

C c_add(C a, C b) {
    C c = c_alloc(a.nt + b.nt);
    memcpy(c.t, a.t, (size_t)a.nt * sizeof(CT));
    memcpy(c.t + a.nt, b.t, (size_t)b.nt * sizeof(CT));
    return c_canon(c);
}

C c_scale(C a, Q k) {
    if (q_sign(k) == 0) return c_zero();
    C c = c_alloc(a.nt);
    for (int i = 0; i < a.nt; i++) { c.t[i] = a.t[i]; c.t[i].k = q_mul(a.t[i].k, k); }
    return c;
}

C c_neg(C a) { return c_scale(a, qi(-1)); }
C c_sub(C a, C b) { return c_add(a, c_neg(b)); }

C c_mul(C a, C b) {
    C c = c_alloc(a.nt * b.nt);
    int w = 0;
    for (int i = 0; i < a.nt; i++)
        for (int j = 0; j < b.nt; j++) {
            CT t;
            t.k = q_mul(a.t[i].k, b.t[j].k);
            for (int l = 0; l < NM_MAXL; l++) t.e[l] = q_add(a.t[i].e[l], b.t[j].e[l]);
            c.t[w++] = t;
        }
    c.nt = w;
    return c_canon(c);
}

int c_is_monomial(C a) { return a.nt == 1; }

/* a / b, when b is a single term */
C c_div(C a, C b) {
    if (b.nt == 0) nm_fail("division by zero");
    if (b.nt > 1) nm_fail("division by %s, which is not a single term (this comes later)", c_to_str(b));
    CT inv;
    inv.k = q_div(qi(1), b.t[0].k);
    for (int l = 0; l < NM_MAXL; l++) inv.e[l] = q_neg(b.t[0].e[l]);
    C m = c_alloc(1);
    m.t[0] = inv;
    return c_mul(a, m);
}

C c_pow_int(C a, int64_t e) {
    if (e < 0) return c_div(c_const(qi(1)), c_pow_int(a, -e));
    C r = c_const(qi(1));
    while (e) {
        if (e & 1) r = c_mul(r, a);
        e >>= 1;
        if (e) a = c_mul(a, a);
    }
    return r;
}

/* a^alpha: exact for a whole alpha, or for a single term whose number has a rational alpha-th power */
int c_pow_q(C a, Q alpha, C *out) {
    if (q_is_int(alpha)) {
        int64_t e;
        if (!z_fits_i64(alpha.num, &e) || e > 100000 || e < -100000) nm_fail("exponent too large");
        *out = c_pow_int(a, e);
        return 1;
    }
    if (a.nt != 1) return 0;
    int64_t p, q;
    if (!z_fits_i64(alpha.num, &p) || !z_fits_i64(alpha.den, &q) || q > 1000) nm_fail("exponent too large");
    Q k;
    if (!q_root_exact(q_pow(a.t[0].k, p), q, &k)) return 0;
    C r = c_alloc(1);
    r.t[0].k = k;
    for (int l = 0; l < NM_MAXL; l++) r.t[0].e[l] = q_mul(a.t[0].e[l], alpha);
    *out = r;
    return 1;
}

int c_equal(C a, C b) { return c_is_zero(c_sub(a, b)); }

C c_persist(C a) {
    C r = a;
    r.t = perm_alloc((size_t)(a.nt > 0 ? a.nt : 1) * sizeof(CT));
    for (int i = 0; i < a.nt; i++) {
        r.t[i].k = q_persist(a.t[i].k);
        for (int l = 0; l < NM_MAXL; l++) r.t[i].e[l] = q_persist(a.t[i].e[l]);
    }
    return r;
}

/* the part of a free of letter idx, and the exponents of idx that occur */
C c_coeff_of(C a, int idx, Q e) {
    C c = c_alloc(a.nt);
    int w = 0;
    for (int i = 0; i < a.nt; i++)
        if (q_cmp(a.t[i].e[idx], e) == 0) { c.t[w] = a.t[i]; c.t[w].e[idx] = q0(); w++; }
    c.nt = w;
    return c;
}

/* ---------------- writing ---------------- */

static char *letter_power(int l, Q e) {           /* a, a^2, a^(1/2) */
    char *s = arena_alloc(strlen(letter_names[l]) + 64);
    if (q_is_int(e) && z_is_one(e.num)) sprintf(s, "%s", letter_names[l]);
    else if (q_is_int(e)) sprintf(s, "%s^%s", letter_names[l], q_to_str(e));
    else sprintf(s, "%s^(%s)", letter_names[l], q_to_str(e));
    return s;
}

/* one term written by hand: 131x^3/(512a^2), -x/4, 2ax, a^(1/2)x^(1/2).
 * extra: an extra letter power appended (the series letter), may be NULL. */
char *ct_str(CT t, const char *extra_name, Q extra_e, int first) {
    int neg = q_sign(t.k) < 0;
    Q k = neg ? q_neg(t.k) : t.k;
    size_t cap = 256;
    for (int l = 0; l < nletters; l++) cap += strlen(letter_names[l]) + 32;
    cap += strlen(q_to_str(k)) + (extra_name ? strlen(extra_name) + 32 : 0);
    char *num = arena_alloc(cap), *den = arena_alloc(cap), *out = arena_alloc(3 * cap);
    char *pn = num, *pd = den;
    *num = *den = 0;
    int nfac = 0, dfac = 0;
    int order[NM_MAXL];                           /* letters in alphabetical order: 2ax, not 2xa */
    for (int l = 0; l < nletters; l++) order[l] = l;
    for (int a = 1; a < nletters; a++)
        for (int b = a; b > 0 && strcmp(letter_names[order[b - 1]], letter_names[order[b]]) > 0; b--) {
            int t2 = order[b]; order[b] = order[b - 1]; order[b - 1] = t2;
        }
    for (int oi = 0; oi < nletters; oi++) {
        int l = order[oi];
        int s = q_sign(t.e[l]);
        if (s > 0) { pn += sprintf(pn, "%s", letter_power(l, t.e[l])); nfac++; }
        if (s < 0) { pd += sprintf(pd, "%s", letter_power(l, q_neg(t.e[l]))); dfac++; }
    }
    if (extra_name && q_sign(extra_e) > 0) {
        char *s = arena_alloc(strlen(extra_name) + 64);
        if (q_is_int(extra_e) && z_is_one(extra_e.num)) sprintf(s, "%s", extra_name);
        else if (q_is_int(extra_e)) sprintf(s, "%s^%s", extra_name, q_to_str(extra_e));
        else sprintf(s, "%s^(%s)", extra_name, q_to_str(extra_e));
        pn += sprintf(pn, "%s", s); nfac++;
    }
    if (extra_name && q_sign(extra_e) < 0) {
        Q e = q_neg(extra_e);
        char *s = arena_alloc(strlen(extra_name) + 64);
        if (q_is_int(e) && z_is_one(e.num)) sprintf(s, "%s", extra_name);
        else if (q_is_int(e)) sprintf(s, "%s^%s", extra_name, q_to_str(e));
        else sprintf(s, "%s^(%s)", extra_name, q_to_str(e));
        pd += sprintf(pd, "%s", s); dfac++;
    }
    char *kn = z_to_str(k.num), *kd = z_to_str(k.den);
    int kd_one = z_is_one(k.den), kn_one = z_is_one(k.num);
    char *o = out;
    if (first) o += sprintf(o, neg ? "-" : "");
    else o += sprintf(o, neg ? " - " : " + ");
    if (!nfac || !kn_one) o += sprintf(o, "%s", kn);
    o += sprintf(o, "%s", num);
    int dparts = (kd_one ? 0 : 1) + dfac;
    if (dparts) {
        int paren = dparts > 1 || (dfac == 1 && !kd_one) || (dfac == 1 && strchr(den, '^'));
        o += sprintf(o, "/%s%s%s%s", paren ? "(" : "", kd_one ? "" : kd, den, paren ? ")" : "");
    }
    return out;
}

/* terms in a readable order: higher total degree first */
static int cmp_for_print(const void *x, const void *y) {
    const CT *a = x, *b = y;
    Q da = q0(), db = q0();
    for (int l = 0; l < NM_MAXL; l++) { da = q_add(da, a->e[l]); db = q_add(db, b->e[l]); }
    int c = q_cmp(db, da);
    if (c) return c;
    for (int l = 0; l < NM_MAXL; l++) {
        c = q_cmp(b->e[l], a->e[l]);
        if (c) return c;
    }
    return 0;
}

char *c_to_str(C a) {
    if (a.nt <= 0) { char *z = arena_alloc(2); strcpy(z, "0"); return z; }
    C s = c_alloc(a.nt);
    memcpy(s.t, a.t, (size_t)a.nt * sizeof(CT));
    qsort(s.t, (size_t)s.nt, sizeof(CT), cmp_for_print);
    size_t cap = 16;
    char **parts = arena_alloc((size_t)s.nt * sizeof(char *));
    for (int i = 0; i < s.nt; i++) { parts[i] = ct_str(s.t[i], NULL, q0(), i == 0); cap += strlen(parts[i]); }
    char *out = arena_alloc(cap), *o = out;
    for (int i = 0; i < s.nt; i++) o += sprintf(o, "%s", parts[i]);
    return out;
}

/* the same, with the series letter's power attached to every term of a */
char *c_term_str(C a, const char *var, Q e, int first) {
    if (a.nt == 1) return ct_str(a.t[0], var, e, first);
    /* several terms: (a + b)x^2 */
    CT one = term_one();
    char *inner = c_to_str(a);
    char *vp = ct_str(one, var, e, 1);
    char *out = arena_alloc(strlen(inner) + strlen(vp) + 16);
    if (q_sign(e) == 0) sprintf(out, "%s(%s)", first ? "" : " + ", inner);
    else sprintf(out, "%s(%s)%s", first ? "" : " + ", inner, vp);
    return out;
}
