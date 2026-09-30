/* Integers in blocks of nine decimal digits, as Newton pointed a number into groups before extracting its root.
 * Values are immutable: every operation returns a fresh number, so inputs may be shared freely. */
#include "nm.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------- memory ---------------- */

typedef struct Chunk { struct Chunk *next; size_t used, size; char data[]; } Chunk;
static Chunk *arena_head;

void *arena_alloc(size_t bytes) {
    bytes = (bytes + 15) & ~(size_t)15;
    if (!arena_head || arena_head->used + bytes > arena_head->size) {
        size_t size = bytes > (1u << 20) ? bytes : (1u << 20);
        Chunk *c = malloc(sizeof(Chunk) + size);
        if (!c) { fprintf(stderr, "out of memory\n"); exit(2); }
        c->next = arena_head; c->used = 0; c->size = size;
        arena_head = c;
    }
    void *p = arena_head->data + arena_head->used;
    arena_head->used += bytes;
    return p;
}

void arena_reset(void) {
    while (arena_head && arena_head->next) { Chunk *n = arena_head->next; free(arena_head); arena_head = n; }
    if (arena_head) arena_head->used = 0;
}

void *perm_alloc(size_t bytes) {
    void *p = malloc(bytes ? bytes : 1);
    if (!p) { fprintf(stderr, "out of memory\n"); exit(2); }
    return p;
}

/* ---------------- errors ---------------- */

jmp_buf nm_on_error;
char nm_error_msg[512];

void nm_fail(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(nm_error_msg, sizeof nm_error_msg, fmt, ap);
    va_end(ap);
    longjmp(nm_on_error, 1);
}

/* ---------------- basics ---------------- */

static Z z_alloc(int n) {
    Z r;
    r.s = 1; r.n = n;
    r.d = n ? arena_alloc((size_t)n * sizeof(uint32_t)) : NULL;
    if (n) memset(r.d, 0, (size_t)n * sizeof(uint32_t));
    return r;
}

static Z z_trim(Z a) {
    while (a.n > 0 && a.d[a.n - 1] == 0) a.n--;
    if (a.n == 0) a.s = 0;
    return a;
}

Z z_zero(void) { Z r = {0, 0, NULL}; return r; }

Z z_from_i64(int64_t v) {
    uint64_t u = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
    Z r = z_alloc(3);
    for (int i = 0; i < 3; i++) { r.d[i] = (uint32_t)(u % NM_BASE); u /= NM_BASE; }
    r.s = v < 0 ? -1 : 1;
    return z_trim(r);
}

Z z_from_dec(const char *s, size_t len) {
    while (len > 1 && *s == '0') { s++; len--; }
    int n = (int)((len + NM_BASE_DIGITS - 1) / NM_BASE_DIGITS);
    Z r = z_alloc(n);
    for (int i = 0; i < n; i++) {
        size_t end = len - (size_t)i * NM_BASE_DIGITS;
        size_t start = end >= NM_BASE_DIGITS ? end - NM_BASE_DIGITS : 0;
        uint32_t v = 0;
        for (size_t k = start; k < end; k++) v = v * 10 + (uint32_t)(s[k] - '0');
        r.d[i] = v;
    }
    return z_trim(r);
}

char *z_to_str(Z a) {
    if (a.s == 0) { char *z = arena_alloc(2); strcpy(z, "0"); return z; }
    char *out = arena_alloc((size_t)a.n * NM_BASE_DIGITS + 2);
    char *p = out;
    if (a.s < 0) *p++ = '-';
    p += sprintf(p, "%u", a.d[a.n - 1]);
    for (int i = a.n - 2; i >= 0; i--) p += sprintf(p, "%09u", a.d[i]);
    return out;
}

int z_is_zero(Z a) { return a.s == 0; }
int z_is_one(Z a) { return a.s == 1 && a.n == 1 && a.d[0] == 1; }

int z_fits_i64(Z a, int64_t *out) {
    if (a.n > 3) return 0;
    uint64_t v = 0;
    for (int i = a.n - 1; i >= 0; i--) {
        if (v > (UINT64_MAX - a.d[i]) / NM_BASE) return 0;
        v = v * NM_BASE + a.d[i];
    }
    if (v > (uint64_t)INT64_MAX) return 0;
    *out = a.s < 0 ? -(int64_t)v : (int64_t)v;
    return 1;
}

int z_cmpabs(Z a, Z b) {
    if (a.n != b.n) return a.n < b.n ? -1 : 1;
    for (int i = a.n - 1; i >= 0; i--)
        if (a.d[i] != b.d[i]) return a.d[i] < b.d[i] ? -1 : 1;
    return 0;
}

int z_cmp(Z a, Z b) {
    if (a.s != b.s) return a.s < b.s ? -1 : 1;
    int c = z_cmpabs(a, b);
    return a.s >= 0 ? c : -c;
}

Z z_neg(Z a) { a.s = -a.s; return a; }
Z z_abs(Z a) { if (a.s < 0) a.s = 1; return a; }

int64_t z_digits(Z a) {
    if (a.s == 0) return 0;
    int64_t k = (int64_t)(a.n - 1) * NM_BASE_DIGITS;
    for (uint32_t t = a.d[a.n - 1]; t; t /= 10) k++;
    return k;
}

Z z_persist(Z a) {
    Z r = a;
    if (a.n) {
        r.d = perm_alloc((size_t)a.n * sizeof(uint32_t));
        memcpy(r.d, a.d, (size_t)a.n * sizeof(uint32_t));
    }
    return r;
}

/* ---------------- addition ---------------- */

static Z add_abs(Z a, Z b) {
    if (a.n < b.n) { Z t = a; a = b; b = t; }
    Z r = z_alloc(a.n + 1);
    uint32_t carry = 0;
    for (int i = 0; i < a.n; i++) {
        uint32_t s = a.d[i] + (i < b.n ? b.d[i] : 0) + carry;
        carry = s >= NM_BASE;
        r.d[i] = carry ? s - NM_BASE : s;
    }
    r.d[a.n] = carry;
    return z_trim(r);
}

static Z sub_abs(Z a, Z b) {                   /* |a| >= |b| */
    Z r = z_alloc(a.n);
    int64_t borrow = 0;
    for (int i = 0; i < a.n; i++) {
        int64_t t = (int64_t)a.d[i] - (i < b.n ? b.d[i] : 0) - borrow;
        borrow = t < 0;
        r.d[i] = (uint32_t)(borrow ? t + NM_BASE : t);
    }
    return z_trim(r);
}

Z z_add(Z a, Z b) {
    if (a.s == 0) return b;
    if (b.s == 0) return a;
    Z r;
    if (a.s == b.s) { r = add_abs(a, b); r.s = a.s; return r; }
    int c = z_cmpabs(a, b);
    if (c == 0) return z_zero();
    if (c > 0) { r = sub_abs(a, b); r.s = a.s; }
    else { r = sub_abs(b, a); r.s = b.s; }
    return r;
}

Z z_sub(Z a, Z b) { return z_add(a, z_neg(b)); }

/* ---------------- multiplication ---------------- */

Z z_mul_small(Z a, uint32_t m) {
    if (a.s == 0 || m == 0) return z_zero();
    Z r = z_alloc(a.n + 1);
    uint64_t carry = 0;
    for (int i = 0; i < a.n; i++) {
        uint64_t t = (uint64_t)a.d[i] * m + carry;
        r.d[i] = (uint32_t)(t % NM_BASE);
        carry = t / NM_BASE;
    }
    r.d[a.n] = (uint32_t)carry;
    r.s = a.s;
    return z_trim(r);
}

Z z_mul(Z a, Z b) {
    if (a.s == 0 || b.s == 0) return z_zero();
    Z r = z_alloc(a.n + b.n);
    for (int i = 0; i < a.n; i++) {
        uint64_t carry = 0, ai = a.d[i];
        if (!ai) continue;
        for (int j = 0; j < b.n; j++) {
            uint64_t t = r.d[i + j] + ai * b.d[j] + carry;
            r.d[i + j] = (uint32_t)(t % NM_BASE);
            carry = t / NM_BASE;
        }
        for (int k = i + b.n; carry; k++) {
            uint64_t t = r.d[k] + carry;
            r.d[k] = (uint32_t)(t % NM_BASE);
            carry = t / NM_BASE;
        }
    }
    r.s = a.s * b.s;
    return z_trim(r);
}

Z z_pow(Z a, unsigned e) {
    Z r = z_from_i64(1);
    while (e) {
        if (e & 1) r = z_mul(r, a);
        e >>= 1;
        if (e) a = z_mul(a, a);
    }
    return r;
}

static const uint32_t POW10[10] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};

Z z_pow10(int64_t k) {
    if (k < 0) nm_fail("internal: negative power of ten");
    int blocks = (int)(k / NM_BASE_DIGITS);
    Z r = z_alloc(blocks + 1);
    r.d[blocks] = POW10[k % NM_BASE_DIGITS];
    return r;
}

Z z_mul_pow10(Z a, int64_t k) {
    if (a.s == 0 || k == 0) return a;
    if (k < 0) nm_fail("internal: negative shift");
    int blocks = (int)(k / NM_BASE_DIGITS);
    Z r = z_alloc(a.n + blocks);
    memcpy(r.d + blocks, a.d, (size_t)a.n * sizeof(uint32_t));
    r.s = a.s;
    return z_mul_small(r, POW10[k % NM_BASE_DIGITS]);
}

/* ---------------- division ---------------- */

static Z divmod_small(Z a, uint32_t m, uint32_t *rem) {
    Z q = z_alloc(a.n);
    uint64_t r = 0;
    for (int i = a.n - 1; i >= 0; i--) {
        uint64_t cur = r * NM_BASE + a.d[i];
        q.d[i] = (uint32_t)(cur / m);
        r = cur % m;
    }
    *rem = (uint32_t)r;
    q.s = a.s;
    return z_trim(q);
}

/* Long division on blocks (Knuth's algorithm D), the same pattern as dividing decimals by hand:
 * estimate a quotient block from the leading blocks, subtract, correct by at most two. */
static void divmod_abs(Z a, Z b, Z *q, Z *r) {
    int n = b.n, m = a.n - b.n;
    uint32_t f = NM_BASE / (b.d[n - 1] + 1);            /* normalise so the divisor's top block is large */
    Z u = z_mul_small(z_abs(a), f), v = z_mul_small(z_abs(b), f);
    uint32_t *ud = arena_alloc((size_t)(a.n + 1) * sizeof(uint32_t));
    memset(ud, 0, (size_t)(a.n + 1) * sizeof(uint32_t));
    memcpy(ud, u.d, (size_t)u.n * sizeof(uint32_t));
    const uint32_t *vd = v.d;
    Z qq = z_alloc(m + 1);
    for (int j = m; j >= 0; j--) {
        uint64_t num = (uint64_t)ud[j + n] * NM_BASE + ud[j + n - 1];
        uint64_t qhat = num / vd[n - 1], rhat = num % vd[n - 1];
        while (qhat >= NM_BASE || (n >= 2 && qhat * vd[n - 2] > rhat * NM_BASE + ud[j + n - 2])) {
            qhat--; rhat += vd[n - 1];
            if (rhat >= NM_BASE) break;
        }
        int64_t borrow = 0;
        uint64_t carry = 0;
        for (int i = 0; i < n; i++) {
            uint64_t p = qhat * vd[i] + carry;
            carry = p / NM_BASE;
            int64_t t = (int64_t)ud[i + j] - (int64_t)(p % NM_BASE) - borrow;
            borrow = t < 0;
            ud[i + j] = (uint32_t)(borrow ? t + NM_BASE : t);
        }
        int64_t t = (int64_t)ud[j + n] - (int64_t)carry - borrow;
        borrow = t < 0;
        ud[j + n] = (uint32_t)(borrow ? t + NM_BASE : t);
        if (borrow) {                                   /* the estimate was one too large: add back */
            qhat--;
            uint32_t c = 0;
            for (int i = 0; i < n; i++) {
                uint32_t s = ud[i + j] + vd[i] + c;
                c = s >= NM_BASE;
                ud[i + j] = c ? s - NM_BASE : s;
            }
            ud[j + n] = (ud[j + n] + c) % NM_BASE;
        }
        qq.d[j] = (uint32_t)qhat;
    }
    Z rr = {1, n, ud};
    uint32_t rem;
    *q = z_trim(qq);
    *r = divmod_small(z_trim(rr), f, &rem);
}

void z_divmod(Z a, Z b, Z *q, Z *r) {
    if (b.s == 0) nm_fail("division by zero");
    Z qq, rr;
    if (z_cmpabs(a, b) < 0) { qq = z_zero(); rr = a; }
    else if (b.n == 1) {
        uint32_t rem;
        qq = divmod_small(z_abs(a), b.d[0], &rem);
        rr = z_from_i64(rem);
    } else divmod_abs(a, b, &qq, &rr);
    if (qq.s) qq.s = a.s * b.s;
    if (rr.s) rr.s = a.s;
    /* the direct operation checks the inverse one: q*b + r must give back a */
    if (z_cmp(z_add(z_mul(qq, b), rr), a) != 0 || z_cmpabs(rr, b) >= 0)
        nm_fail("internal check failed in division (vitiose)");
    *q = qq; *r = rr;
}

Z z_div_round(Z a, Z b) {
    Z q, r;
    z_divmod(a, b, &q, &r);
    if (z_cmpabs(z_mul_small(r, 2), b) >= 0) q = z_add(q, z_from_i64(a.s * b.s));
    return q;
}

Z z_gcd(Z a, Z b) {
    a = z_abs(a); b = z_abs(b);
    while (b.s) { Z q, r; z_divmod(a, b, &q, &r); a = b; b = r; }
    return a;
}

/* floor of the r-th root, by Newton's iteration on integers from a starting value above the root */
Z z_iroot(Z a, unsigned r) {
    if (a.s < 0) nm_fail("internal: root of a negative integer");
    if (a.s == 0 || r == 1) return a;
    Z x = z_pow10((z_digits(a) + r - 1) / r);
    for (;;) {
        Z q, rem;
        z_divmod(a, z_pow(x, r - 1), &q, &rem);
        Z y = z_add(z_mul_small(x, r - 1), q);
        z_divmod(y, z_from_i64(r), &y, &rem);
        if (z_cmp(y, x) >= 0) return x;
        x = y;
    }
}
