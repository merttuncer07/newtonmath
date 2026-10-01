/* Quantities in letters: sums of terms  k · a^e1 · b^e2 ...  with rational k and rational exponents, so that
 * a^(1/2) and a^(-1) are written as Newton wrote them, "ob analogiam rei". Division is by a single term only,
 * as in all of Newton's examples (xx/64a, 131x^3/512aa); anything else is refused with the reason. */
#include "nm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------- letters ---------------- */

static char *letter_names[NM_MAXL];        /* how a letter is written */
static char *letter_keys[NM_MAXL];         /* how it is looked up: the same, or "#..." for a surd */
static int nletters;

/* A surd is a letter with its own equation (Methodus, Problem 1, Example 3: "pro singulis pono totidem literas"):
 *   r^d = e[0] + e[1] r + ... + e[d-1] r^(d-1),
 * and a choice of root: a certified bracket, a radical of a positive quantity, or i. */
typedef struct { int d; C *e; Root *root; C radA; int radn; int radneg; int imag;
                 Ball (*value)(void *data, int64_t prec); void *data; } Alg;   /* value: a named number, no equation */
static Alg alg[NM_MAXL];

int letter_index(const char *name, size_t len) {
    for (int i = 0; i < nletters; i++)
        if (strlen(letter_keys[i]) == len && !strncmp(letter_keys[i], name, len)) return i;
    if (nletters == NM_MAXL) nm_fail("too many different letters and surds (at most %d)", NM_MAXL);
    char *s = perm_alloc(len + 1);
    memcpy(s, name, len); s[len] = 0;
    letter_names[nletters] = letter_keys[nletters] = s;
    memset(&alg[nletters], 0, sizeof(Alg));
    return nletters++;
}

const char *letter_name(int i) { return letter_names[i]; }
int letter_is_surd(int i) { return alg[i].d > 0; }
int letter_is_named(int i) { return alg[i].value != NULL; }
Ball named_ball(int i, int64_t prec) { return alg[i].value(alg[i].data, prec); }

/* a number known by its value only, kept as a letter: exp(1), sin(1) ... ("pro singulis pono literas") */
C c_named(const char *disp, Ball (*value)(void *, int64_t), void *data) {
    /* The callback identifies the meaning as well as the display: a user rule
     * called atan must not replace the circular area with the same spelling. */
    for (int l = 0; l < nletters; l++) if (alg[l].value == value && !strcmp(letter_names[l], disp)) return c_letter(l);
    char key[32];
    snprintf(key, sizeof key, "#n%d", nletters);
    int l = letter_index(key, strlen(key));
    size_t n = strlen(disp);
    letter_names[l] = perm_alloc(n + 1);
    memcpy(letter_names[l], disp, n + 1);
    alg[l].value = value; alg[l].data = data;
    return c_letter(l);
}
int letter_is_imag(int i) { return alg[i].imag; }
Root *surd_root(int i) { return alg[i].root; }
int surd_radical(int i, C *A, int *n, int *neg) { *A = alg[i].radA; *n = alg[i].radn; *neg = alg[i].radneg; return alg[i].radn > 0; }
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

static C c_alloc(int nt);
static C c_pow_raw(int l, int k);

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

static C c_pow_raw(int l, int k) {            /* r^k, unreduced */
    C c = c_alloc(1);
    c.t[0] = term_one();
    c.t[0].e[l] = qi(k);
    return c;
}

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

static C c_mul_raw(C a, C b);

/* replace r^k (k >= d) by the equation of r, until every surd is below its degree */
static C c_reduce(C a) {
    for (int guard = 0; guard < 100000; guard++) {
        int hit = -1, l = -1;
        for (int t = 0; t < a.nt && hit < 0; t++)
            for (int m = 0; m < nletters; m++)
                if (alg[m].d && q_cmp(a.t[t].e[m], qi(alg[m].d)) >= 0) { hit = t; l = m; break; }
        if (hit < 0) return a;
        CT t = a.t[hit];
        t.e[l] = q_sub(t.e[l], qi(alg[l].d));
        C mono; mono.nt = 1; mono.t = arena_alloc(sizeof(CT)); mono.t[0] = t;
        C rest = c_alloc(a.nt - 1);
        int w = 0;
        for (int u = 0; u < a.nt; u++) if (u != hit) rest.t[w++] = a.t[u];
        rest.nt = w;
        C sub = c_zero();                          /* e[0] + e[1] r + ... */
        for (int i = 0; i < alg[l].d; i++) sub = c_add(sub, c_mul_raw(alg[l].e[i], c_pow_raw(l, i)));
        a = c_add(rest, c_mul_raw(mono, sub));
    }
    nm_fail("internal: reduction by a surd's equation did not end");
}

C c_mul(C a, C b) { return c_reduce(c_mul_raw(a, b)); }

static C c_mul_raw(C a, C b) {
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

int c_has_plain(C a) {                        /* a letter that is not a surd */
    for (int t = 0; t < a.nt; t++)
        for (int l = 0; l < nletters; l++) if (!alg[l].d && !alg[l].value && q_sign(a.t[t].e[l])) return 1;
    return 0;
}

static int mono_has_surd(CT t) {
    for (int l = 0; l < nletters; l++) if (alg[l].d && q_sign(t.e[l])) return 1;
    return 0;
}

/* ---------------- inverses: clearing a surd from a denominator ---------------- */

C c_inv(C a);

typedef struct { int deg; C *c; } RP;          /* a polynomial in one surd r, coefficients below it */

static RP rp_alloc(int deg) {
    RP p; p.deg = deg; p.c = arena_alloc((size_t)(deg + 2) * sizeof(C));
    for (int i = 0; i <= deg; i++) p.c[i] = c_zero();
    return p;
}
static RP rp_trim(RP p) { while (p.deg >= 0 && c_is_zero(p.c[p.deg])) p.deg--; return p; }
static RP rp_from(C a, int r) {                /* a as a polynomial in r (whole powers) */
    int deg = 0;
    for (int t = 0; t < a.nt; t++) { int64_t e; z_fits_i64(a.t[t].e[r].num, &e); if (e > deg) deg = (int)e; }
    RP p = rp_alloc(deg);
    for (int t = 0; t < a.nt; t++) {
        int64_t e; z_fits_i64(a.t[t].e[r].num, &e);
        CT u = a.t[t]; u.e[r] = qi(0);
        C one; one.nt = 1; one.t = arena_alloc(sizeof(CT)); one.t[0] = u;
        p.c[e] = c_add(p.c[e], one);
    }
    return rp_trim(p);
}
static C rp_to(RP p, int r) {
    C out = c_zero();
    for (int i = 0; i <= p.deg; i++) out = c_add(out, c_mul(p.c[i], c_pow_raw(r, i)));
    return c_reduce(out);
}
static RP rp_sub(RP a, RP b) {
    RP r = rp_alloc(a.deg > b.deg ? a.deg : b.deg);
    for (int i = 0; i <= r.deg; i++) r.c[i] = c_sub(i <= a.deg ? a.c[i] : c_zero(), i <= b.deg ? b.c[i] : c_zero());
    return rp_trim(r);
}
static RP rp_mul(RP a, RP b) {
    if (a.deg < 0 || b.deg < 0) return rp_alloc(-1);
    RP r = rp_alloc(a.deg + b.deg);
    for (int i = 0; i <= a.deg; i++) for (int j = 0; j <= b.deg; j++) r.c[i + j] = c_add(r.c[i + j], c_mul(a.c[i], b.c[j]));
    return rp_trim(r);
}
static void rp_divmod(RP a, RP b, RP *q, RP *rem) {
    RP qq = rp_alloc(a.deg - b.deg >= 0 ? a.deg - b.deg : 0);
    C inv = c_inv(b.c[b.deg]);
    while (a.deg >= b.deg && a.deg >= 0) {
        int sh = a.deg - b.deg;
        C f = c_mul(a.c[a.deg], inv);
        qq.c[sh] = c_add(qq.c[sh], f);
        RP m = rp_alloc(b.deg + sh);
        for (int i = 0; i <= b.deg; i++) m.c[i + sh] = c_mul(b.c[i], f);
        RP na = rp_sub(a, m);
        na.c[a.deg] = c_zero();                        /* exact cancellation of the leading term */
        a = rp_trim(na);
    }
    *q = rp_trim(qq); *rem = a;
}

/* 1/a: a single term directly; otherwise, for the latest surd r in a, Euclid between a and r's equation */
C c_inv(C a) {
    if (a.nt == 0) nm_fail("division by zero");
    if (a.nt == 1 && !mono_has_surd(a.t[0])) {
        CT inv;
        inv.k = q_div(qi(1), a.t[0].k);
        for (int l = 0; l < NM_MAXL; l++) inv.e[l] = q_neg(a.t[0].e[l]);
        C m = c_alloc(1);
        m.t[0] = inv;
        return m;
    }
    /* a common single term in the ordinary letters comes out first */
    CT common = a.t[0];
    for (int l = 0; l < NM_MAXL; l++) if (l >= nletters || alg[l].d) common.e[l] = qi(0);
    for (int t = 1; t < a.nt; t++)
        for (int l = 0; l < nletters; l++)
            if (!alg[l].d && q_cmp(a.t[t].e[l], common.e[l]) != 0)
                nm_fail("division by %s, whose ordinary letters are not a single term (this comes later)", c_to_str(a));
    common.k = qi(1);
    C M; M.nt = 1; M.t = arena_alloc(sizeof(CT)); M.t[0] = common;
    C rest = c_mul_raw(a, c_inv(M));
    if (!mono_has_surd(M.t[0]) && M.nt == 1) {
        int r = -1;
        for (int t = 0; t < rest.nt; t++) for (int l = 0; l < nletters; l++) if (alg[l].d && q_sign(rest.t[t].e[l]) && l > r) r = l;
        if (r < 0) return c_mul(c_inv(rest), c_inv(M));
        RP A = rp_alloc(alg[r].d);                       /* r^d - e[d-1] r^(d-1) - ... - e[0] */
        A.c[alg[r].d] = c_const(qi(1));
        for (int i = 0; i < alg[r].d; i++) A.c[i] = c_neg(alg[r].e[i]);
        RP r0 = A, r1 = rp_from(rest, r), t0 = rp_alloc(-1), t1 = rp_alloc(0);
        t1.c[0] = c_const(qi(1)); t1 = rp_trim(t1);
        while (r1.deg > 0) {
            RP q, rm;
            rp_divmod(r0, r1, &q, &rm);
            RP t2 = rp_sub(t0, rp_mul(q, t1));
            r0 = r1; r1 = rm; t0 = t1; t1 = t2;
        }
        if (r1.deg < 0) {
            /* the equation of r factors: keep the factor its root satisfies, and divide again */
            Root *rt = alg[r].root;
            int rational = 1;
            for (int i = 0; i <= r0.deg; i++) { Q k; if (!c_const_value(r0.c[i], &k)) rational = 0; }
            if (!rt || !rational)
                nm_fail("the equation of %s factors, and this denominator vanishes with one factor (define the number by the factor it satisfies)",
                        letter_names[r]);
            RP h, rm;
            rp_divmod(A, r0, &h, &rm);
            RP g = r0;
            Q lo = q_make(z_sub(rt->X, z_from_i64(rt->w)), z_pow10(rt->D)), hi = q_make(z_add(rt->X, z_from_i64(rt->w)), z_pow10(rt->D));
            Q glo = qi(0), ghi = qi(0), k;
            for (int i = g.deg; i >= 0; i--) {
                c_const_value(g.c[i], &k);
                glo = q_add(q_mul(glo, lo), k); ghi = q_add(q_mul(ghi, hi), k);
            }
            int in_g = rt->w == 0 ? q_sign(glo) == 0 : q_sign(glo) * q_sign(ghi) <= 0;
            RP keep = in_g ? g : h;
            C lead = c_inv(keep.c[keep.deg]);
            alg[r].d = keep.deg;
            alg[r].e = perm_alloc((size_t)keep.deg * sizeof(C));
            for (int i = 0; i < keep.deg; i++) alg[r].e[i] = c_persist(c_neg(c_mul(keep.c[i], lead)));
            return c_inv(c_reduce(a));
        }
        C inv = c_mul(rp_to(t1, r), c_inv(r1.c[0]));
        return c_mul(inv, c_inv(M));
    }
    return c_mul(c_inv(rest), c_inv(M));
}

/* a / b: by a single term directly, by anything else through its inverse */
C c_div(C a, C b) {
    if (b.nt == 0) nm_fail("division by zero");
    return c_mul(a, c_inv(b));
}

C c_pow_int(C a, int64_t e) {
    if (e < 0) return c_inv(c_pow_int(a, -e));
    C r = c_const(qi(1));
    while (e) {
        if (e & 1) r = c_mul(r, a);
        e >>= 1;
        if (e) a = c_mul(a, a);
    }
    return r;
}

/* ---------------- surds: new letters with their equations ---------------- */

static int surd_find(const char *disp) {
    for (int l = 0; l < nletters; l++) if (alg[l].d && !strcmp(letter_names[l], disp)) return l;
    return -1;
}

static int surd_make(const char *disp, int d, C *e);
static int surd_new(const char *disp, int d, C *e) {
    int f = surd_find(disp);
    if (f >= 0) return f;
    return surd_make(disp, d, e);
}

static int surd_make(const char *disp, int d, C *e) {
    char key[32];
    snprintf(key, sizeof key, "#s%d", nletters);
    int l = letter_index(key, strlen(key));
    size_t n = strlen(disp);
    letter_names[l] = perm_alloc(n + 1);
    memcpy(letter_names[l], disp, n + 1);
    alg[l].d = d;
    alg[l].e = perm_alloc((size_t)d * sizeof(C));
    for (int i = 0; i < d; i++) alg[l].e[i] = c_persist(e[i]);
    return l;
}

/* a root found by Newton's resolution, kept as a letter with its equation: let r = root of ... near ... */
C c_surd_from_root(const char *name, Root *r) {
    root_refine(r, 20);
    if (r->certified && r->w == 0) return c_const(q_make(r->X, z_pow10(r->D)));   /* a rational root */
    C *e = arena_alloc((size_t)r->deg * sizeof(C));
    for (int i = 0; i < r->deg; i++) e[i] = c_const(q_neg(q_make(r->c[i], r->c[r->deg])));
    int l = surd_make(name, r->deg, e);
    alg[l].root = r;
    return c_letter(l);
}

C c_imag_unit(void) {                           /* i, with i^2 = -1 */
    C e[2] = {c_const(qi(-1)), c_zero()};
    int l = surd_new("i", 2, e);
    alg[l].imag = 1;
    return c_letter(l);
}

/* c^(1/n) for a rational c: exact when it is; otherwise k times a surd, with the n-th powers taken out first */
C c_radical_q(Q c, int64_t n) {
    Q ex;
    if (q_root_exact(c, n, &ex)) return c_const(ex);
    int neg = q_sign(c) < 0;
    if (neg && n % 2 == 0) {
        if (n != 2) nm_fail("an even root of a negative number other than the square root comes later");
        return c_mul(c_imag_unit(), c_radical_q(q_neg(c), 2));
    }
    Z A = z_mul(z_abs(c.num), z_pow(c.den, (unsigned)(n - 1)));    /* (A/B)^(1/n) = (A B^(n-1))^(1/n) / B */
    Z out = z_from_i64(1), rem;
    for (uint32_t p = 2; p < 2000; p++) {                         /* take out small n-th powers */
        Z pn = z_pow(z_from_i64(p), (unsigned)n);
        if (z_cmp(pn, A) > 0) break;
        for (;;) {
            Z q;
            z_divmod(A, pn, &q, &rem);
            if (rem.s) break;
            A = q; out = z_mul_small(out, p);
        }
    }
    Q coef = q_make(out, c.den);
    if (neg) coef = q_neg(coef);
    if (z_is_one(A)) return c_const(coef);
    char *num = z_to_str(A);
    char *disp = arena_alloc(strlen(num) + 32);
    if (n == 2) sprintf(disp, "sqrt(%s)", num); else sprintf(disp, "%s^(1/%lld)", num, (long long)n);
    int f = surd_find(disp);
    if (f < 0) {
        C *e = arena_alloc((size_t)n * sizeof(C));
        for (int i = 0; i < n; i++) e[i] = c_zero();
        e[0] = c_const(q_from_z(A));
        f = surd_new(disp, (int)n, e);
        /* its root: the positive real one, by Newton's resolution from the integer root */
        Poly p; p.deg = (int)n; p.var = "y";
        p.c = arena_alloc((size_t)(n + 1) * sizeof(Q));
        for (int i = 0; i <= n; i++) p.c[i] = qi(0);
        p.c[n] = qi(1); p.c[0] = q_neg(q_from_z(A));
        Z g = z_iroot(z_mul_pow10(A, 6 * n), (unsigned)n);
        alg[f].root = root_new(p, q_make(g, z_pow10(6)), 6);
    }
    return c_scale(c_letter(f), coef);
}

/* A^(1/n) for a quantity made of numbers and surds only: a new surd over the earlier ones (a tower) */
C c_radical_c(C A, int64_t n) {
    Q k;
    if (c_const_value(A, &k)) return c_radical_q(k, n);
    if (c_has_plain(A)) nm_fail("internal: a radical of letters");
    if (n > 12) nm_fail("root index too large for a surd over surds");
    char *inner = c_to_str(A);
    char *disp = arena_alloc(strlen(inner) + 32);
    if (n == 2) sprintf(disp, "sqrt(%s)", inner); else sprintf(disp, "(%s)^(1/%lld)", inner, (long long)n);
    int f = surd_find(disp);
    if (f < 0) {
        C *e = arena_alloc((size_t)n * sizeof(C));
        for (int i = 0; i < n; i++) e[i] = c_zero();
        e[0] = A;
        f = surd_new(disp, (int)n, e);
        alg[f].radA = c_persist(A); alg[f].radn = (int)n;
    }
    return c_letter(f);
}

/* a^alpha: exact for a whole alpha, for a single term, or for numbers and surds (then a new surd) */
int c_pow_q(C a, Q alpha, C *out) {
    if (q_is_int(alpha)) {
        int64_t e;
        if (!z_fits_i64(alpha.num, &e) || e > 100000 || e < -100000) nm_fail("exponent too large");
        *out = c_pow_int(a, e);
        return 1;
    }
    int64_t p, q;
    if (!z_fits_i64(alpha.num, &p) || !z_fits_i64(alpha.den, &q) || q > 1000) nm_fail("exponent too large");
    if (a.nt == 1) {
        CT t = a.t[0];
        C number = c_radical_q(q_pow(t.k, p), q);                 /* k^alpha */
        CT plain = t, surd = t;
        plain.k = qi(1); surd.k = qi(1);
        for (int l = 0; l < NM_MAXL; l++) {
            if (l < nletters && alg[l].d) plain.e[l] = qi(0); else surd.e[l] = qi(0);
            plain.e[l] = q_mul(plain.e[l], alpha);
        }
        C P; P.nt = 1; P.t = arena_alloc(sizeof(CT)); P.t[0] = plain;
        C S; S.nt = 1; S.t = arena_alloc(sizeof(CT)); S.t[0] = surd;
        C sp = mono_has_surd(surd) ? c_radical_c(c_pow_int(S, p), q) : c_const(qi(1));
        *out = c_mul(c_mul(number, P), sp);
        return 1;
    }
    if (c_has_plain(a)) return 0;
    *out = c_radical_c(c_pow_int(a, p), q);
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

static char *letter_power(int l, Q e) {           /* a, a^2, a^(1/2), sqrt(2), (2^(1/3))^2 */
    char *s = arena_alloc(strlen(letter_names[l]) + 64);
    if (q_is_int(e) && z_is_one(e.num)) sprintf(s, "%s", letter_names[l]);
    else if (alg[l].d && strspn(letter_names[l], "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789") != strlen(letter_names[l]))
        sprintf(s, "(%s)^%s", letter_names[l], q_to_str(e));
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
    int order[NM_MAXL];                           /* surds first, then letters in alphabetical order: 2sqrt(2)ax */
    for (int l = 0; l < nletters; l++) order[l] = l;
    for (int a = 1; a < nletters; a++)
        for (int b = a; b > 0; b--) {
#define NUMLIKE(l) (alg[l].d > 0 || alg[l].value)
            int sb = NUMLIKE(order[b]) + (NUMLIKE(order[b]) && !alg[order[b]].imag);   /* 2: surds and named, 1: i, 0: letters */
            int sa = NUMLIKE(order[b - 1]) + (NUMLIKE(order[b - 1]) && !alg[order[b - 1]].imag);
            int before = sb > sa || (sb == sa && strcmp(letter_names[order[b - 1]], letter_names[order[b]]) > 0);
            if (!before) break;
            int t2 = order[b]; order[b] = order[b - 1]; order[b - 1] = t2;
        }
    for (int oi = 0; oi < nletters; oi++) {
        int l = order[oi];
        if (alg[l].imag) continue;
        int s = q_sign(t.e[l]);
        if (s > 0) { pn += sprintf(pn, "%s", letter_power(l, t.e[l])); nfac++; }
        if (s < 0) { pd += sprintf(pd, "%s", letter_power(l, q_neg(t.e[l]))); dfac++; }
    }
    for (int oi = 0; oi < nletters; oi++) {          /* i last: sqrt(3)i, 2a*i */
        int l = order[oi];
        if (!alg[l].imag || q_sign(t.e[l]) == 0) continue;
        if (pn > num && ((pn[-1] >= 'a' && pn[-1] <= 'z') || (pn[-1] >= 'A' && pn[-1] <= 'Z'))) *pn++ = '*';
        pn += sprintf(pn, "i"); nfac++;
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

/* Terms: descending total degree, then descending exponents in alphabetical letter order. */
static int print_order[NM_MAXL];

static int cmp_for_print(const void *x, const void *y) {
    const CT *a = x, *b = y;
    Q da = q0(), db = q0();
    for (int l = 0; l < nletters; l++) { da = q_add(da, a->e[l]); db = q_add(db, b->e[l]); }
    int c = q_cmp(db, da);
    if (c) return c;
    for (int i = 0; i < nletters; i++) {
        int l = print_order[i];
        c = q_cmp(b->e[l], a->e[l]);
        if (c) return c;
    }
    return 0;
}

static void set_print_order(void) {
    for (int l = 0; l < nletters; l++) {
        int i = l;
        while (i > 0 && strcmp(letter_names[l], letter_names[print_order[i - 1]]) < 0) {
            print_order[i] = print_order[i - 1];
            i--;
        }
        print_order[i] = l;
    }
}

int c_leading_index(C a) {
    set_print_order();
    int best = -1;
    for (int i = 0; i < a.nt; i++)
        if (best < 0 || cmp_for_print(&a.t[i], &a.t[best]) < 0) best = i;
    return best;
}

char *c_to_str(C a) {
    if (a.nt <= 0) { char *z = arena_alloc(2); strcpy(z, "0"); return z; }
    C s = c_alloc(a.nt);
    memcpy(s.t, a.t, (size_t)a.nt * sizeof(CT));
    set_print_order();
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
