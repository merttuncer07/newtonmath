/* Roots by Newton's resolution, and balls (decimal midpoint with a guaranteed radius).
 *
 * A root keeps its equation and its current approximation. Asking for more places continues from where it
 * stopped, as Newton's residual equation lets the work go on "quousque placuerit". Each pass runs at the
 * precision it can deliver, and the precision doubles from pass to pass. The bracket is confirmed by the direct
 * operation: substitute both ends into the equation exactly and see the sign change. */
#include "nm.h"

#include <string.h>

/* ---------------- roots ---------------- */

/* p(X / 10^D) * 10^(D*deg), exactly, by nested multiplication: (((c_n X + c_{n-1} 10^D) X + ...) */
static Z eval_scaled(const Z *c, int deg, Z X, int64_t D) {
    Z h = c[deg];
    for (int i = deg - 1; i >= 0; i--) h = z_add(z_mul(h, X), z_mul_pow10(c[i], D * (deg - i)));
    return h;
}

static void newton_at(Root *r, Z *dc, int64_t D) {
    for (int it = 0; it < 200; it++) {
        Z h = eval_scaled(r->c, r->deg, r->X, D);
        if (h.s == 0) return;
        Z hp = eval_scaled(dc, r->deg - 1, r->X, D);
        if (hp.s == 0) nm_fail("the derivative vanishes at the approximation; choose another starting value");
        Z delta = z_div_round(h, hp);      /* the correction p(x)/p'(x), in units of 10^-D */
        r->X = z_sub(r->X, delta);
        if (z_cmpabs(delta, z_from_i64(1)) <= 0) return;
    }
    nm_fail("Newton's resolution does not settle near the given value");
}

Root *root_new(Poly p, Q guess, int64_t start_places) {
    if (p.deg < 1 || !p.var) nm_fail("the equation has no unknown");
    Z l = z_from_i64(1);                     /* clear denominators: multiply by their least common multiple */
    for (int i = 0; i <= p.deg; i++) {
        Z g = z_gcd(l, p.c[i].den), q, rem;
        z_divmod(z_mul(l, p.c[i].den), g, &q, &rem);
        l = q;
    }
    Root *r = perm_alloc(sizeof *r);
    r->deg = p.deg;
    r->c = perm_alloc((size_t)(p.deg + 1) * sizeof(Z));
    for (int i = 0; i <= p.deg; i++) {
        Z q, rem;
        z_divmod(z_mul(p.c[i].num, l), p.c[i].den, &q, &rem);
        r->c[i] = z_persist(q);
    }
    size_t vl = strlen(p.var);
    char *v = perm_alloc(vl + 1);
    memcpy(v, p.var, vl + 1);
    r->var = v;
    r->D = start_places < 4 ? 4 : start_places;
    r->X = z_persist(z_div_round(z_mul_pow10(guess.num, r->D), guess.den));
    r->w = 0;
    r->certified = 0;
    return r;
}

static int sign_at(Root *r, Z X, int64_t D) { return eval_scaled(r->c, r->deg, X, D).s; }

void root_refine(Root *r, int64_t places) {
    int64_t target = places + 3;
    if (r->certified && r->D >= target) return;
    Z *dc = arena_alloc((size_t)r->deg * sizeof(Z));
    for (int i = 0; i < r->deg; i++) dc[i] = z_mul_small(r->c[i + 1], (uint32_t)(i + 1));
    int64_t D = r->D;
    newton_at(r, dc, D);
    while (D < target) {                      /* double the places each pass, as far as needed */
        int64_t D2 = 2 * D < target ? 2 * D : target;
        r->X = z_mul_pow10(r->X, D2 - D);
        D = D2;
        newton_at(r, dc, D);
    }
    r->D = D;
    /* the direct operation: the equation must change sign across the bracket */
    r->certified = 0;
    if (sign_at(r, r->X, D) == 0) { r->w = 0; r->certified = 1; }
    else {
        for (int64_t w = 1; w <= 500; w *= 2) {
            Z one = z_from_i64(w);
            if (sign_at(r, z_sub(r->X, one), D) * sign_at(r, z_add(r->X, one), D) < 0) {
                r->w = w; r->certified = 1;
                break;
            }
        }
    }
    r->X = z_persist(r->X);
}

char *root_equation_str(Root *r) {
    Poly p;
    p.deg = r->deg; p.var = r->var;
    p.c = arena_alloc((size_t)(r->deg + 1) * sizeof(Q));
    for (int i = 0; i <= r->deg; i++) p.c[i] = q_from_z(r->c[i]);
    return p_to_str(p);
}

/* ---------------- balls ---------------- */

static Z ceil_div_pow10(Z a, int64_t k) {    /* a >= 0 */
    Z q, rem;
    z_divmod(a, z_pow10(k), &q, &rem);
    return rem.s ? z_add(q, z_from_i64(1)) : q;
}

/* keep `prec` significant digits of the midpoint; the radius absorbs the rounding */
static Ball b_trim(Ball a, int64_t prec) {
    int64_t k = z_digits(a.m) - prec;
    if (k > 0) {
        a.m = z_div_round(a.m, z_pow10(k));
        a.r = z_add(ceil_div_pow10(a.r, k), z_from_i64(1));
        a.e += k;
    }
    return a;
}

Ball b_from_q(Q a, int64_t prec) {
    int64_t s = prec - (z_digits(a.num) - z_digits(a.den)) + 1;
    Ball b;
    Z q, rem;
    if (s >= 0) z_divmod(z_mul_pow10(a.num, s), a.den, &q, &rem);
    else z_divmod(a.num, z_mul(a.den, z_pow10(-s)), &q, &rem);
    b.m = q; b.e = -s; b.r = z_from_i64(rem.s ? 1 : 0);
    return b_trim(b, prec);
}

Ball b_from_root(Root *r, int64_t prec) {
    int64_t mag = z_digits(r->X) - r->D;     /* digits before the decimal point, roughly */
    int64_t places = prec - (mag > 0 ? mag : 0);
    root_refine(r, places);
    if (!r->certified) nm_fail("the root of %s could not be certified (no sign change; a double root?)",
                               root_equation_str(r));
    Ball b = {r->X, -r->D, z_from_i64(r->w)};
    return b_trim(b, prec);
}

static Ball b_align(Ball a, int64_t e) {     /* rewrite with a smaller exponent e <= a.e */
    a.m = z_mul_pow10(a.m, a.e - e);
    a.r = z_mul_pow10(a.r, a.e - e);
    a.e = e;
    return a;
}

static int64_t top_exp(Ball a) {             /* |midpoint| + radius < 10^top_exp */
    int64_t dm = z_digits(a.m), dr = z_digits(a.r);
    return a.e + (dm > dr ? dm : dr) + 1;
}

Ball b_add(Ball a, Ball b, int64_t prec) {
    if (a.e < b.e) { Ball t = a; a = b; b = t; }
    /* b far below a's last digit: it only widens the radius */
    if (a.m.s && a.e - b.e > prec + 20 && top_exp(b) <= a.e) {
        a.r = z_add(a.r, z_from_i64(1));
        return a;
    }
    a = b_align(a, b.e);
    Ball r = {z_add(a.m, b.m), b.e, z_add(a.r, b.r)};
    return b_trim(r, prec);
}

Ball b_sub(Ball a, Ball b, int64_t prec) { b.m = z_neg(b.m); return b_add(a, b, prec); }

Ball b_mul(Ball a, Ball b, int64_t prec) {
    Ball r;
    r.m = z_mul(a.m, b.m);
    r.r = z_add(z_add(z_mul(z_abs(a.m), b.r), z_mul(z_abs(b.m), a.r)), z_mul(a.r, b.r));
    r.e = a.e + b.e;
    return b_trim(r, prec);
}

Ball b_div(Ball a, Ball b, int64_t prec) {
    Z bm = z_abs(b.m);
    if (z_cmp(bm, b.r) <= 0) nm_fail("division by a number that may be zero");
    int64_t s = prec + z_digits(b.m) - z_digits(a.m) + 2;
    if (s < 0) s = 0;
    Z q, rem;
    z_divmod(z_mul_pow10(a.m, s), b.m, &q, &rem);
    /* |x/y - ma/mb| <= (ra |mb| + |ma| rb) / (|mb| (|mb| - rb)) */
    Z num = z_mul_pow10(z_add(z_mul(a.r, bm), z_mul(z_abs(a.m), b.r)), s);
    Z den = z_mul(bm, z_sub(bm, b.r));
    Z rq, rr;
    z_divmod(num, den, &rq, &rr);
    if (rr.s) rq = z_add(rq, z_from_i64(1));
    Ball r = {q, a.e - b.e - s, z_add(rq, z_from_i64(1))};
    return b_trim(r, prec);
}

Ball b_pow(Ball a, int64_t e, int64_t prec) {
    if (e < 0) {
        Ball one = {z_from_i64(1), 0, z_zero()};
        return b_div(one, b_pow(a, -e, prec), prec);
    }
    Ball r = {z_from_i64(1), 0, z_zero()};
    while (e) {
        if (e & 1) r = b_mul(r, a, prec);
        e >>= 1;
        if (e) a = b_mul(a, a, prec);
    }
    return r;
}

/* the largest k <= want such that radius * 10^e <= 0.5 * 10^-k; -1 if not even k = 0 */
int64_t b_guaranteed_places(Ball a, int64_t want) {
    if (a.r.s == 0) return want;
    Z r2 = z_mul_small(a.r, 2);
    int64_t t = z_digits(r2);                /* smallest t with 10^t >= 2r, unless 2r is a power of ten */
    if (z_cmp(z_pow10(t - 1), r2) == 0) t--;
    int64_t k = -a.e - t;
    if (k > want) k = want;
    return k < 0 ? -1 : k;
}

Z b_round_to_places(Ball a, int64_t places) {
    if (a.e + places >= 0) return z_mul_pow10(a.m, a.e + places);
    return z_div_round(a.m, z_pow10(-(a.e + places)));
}

char *fixed_str(Z scaled, int64_t places) {
    char *digits = z_to_str(z_abs(scaled));
    size_t n = strlen(digits);
    size_t need = (size_t)places + 1;
    char *out = arena_alloc((n > need ? n : need) + (size_t)places + 4);
    char *o = out;
    if (scaled.s < 0) *o++ = '-';
    char *padded = arena_alloc((n > need ? n : need) + 1);
    size_t pad = n < need ? need - n : 0;
    memset(padded, '0', pad);
    memcpy(padded + pad, digits, n + 1);
    size_t len = pad + n, intlen = len - (size_t)places;
    memcpy(o, padded, intlen); o += intlen;
    if (places > 0) { *o++ = '.'; memcpy(o, padded + intlen, (size_t)places); o += places; }
    *o = 0;
    return out;
}
