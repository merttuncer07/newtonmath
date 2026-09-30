/* newtonmath: shared declarations.
 *
 * Numbers are kept exact wherever possible. An integer is a list of blocks of nine decimal digits (base 10^9),
 * least significant block first. A rational is a reduced pair of integers. A root is kept as its equation together
 * with a decimal approximation whose bracket has been verified by a sign change. An approximate value is a ball:
 * a decimal midpoint with a guaranteed radius. */
#ifndef NM_H
#define NM_H

#include <stdint.h>
#include <stddef.h>

/* ---- errors and memory ---- */
__attribute__((noreturn, format(printf, 1, 2))) void nm_fail(const char *fmt, ...);            /* abort the current statement with a message */
void *arena_alloc(size_t bytes);               /* temporary memory, freed after each statement */
void arena_reset(void);
void *perm_alloc(size_t bytes);                /* memory that outlives the statement */

/* ---- integers ---- */
#define NM_BASE 1000000000u
#define NM_BASE_DIGITS 9

typedef struct {
    int s;          /* sign: -1, 0, +1 */
    int n;          /* number of blocks in use */
    uint32_t *d;    /* blocks, least significant first, each < NM_BASE */
} Z;

Z z_zero(void);
Z z_from_i64(int64_t v);
Z z_from_dec(const char *digits, size_t len);  /* non-negative decimal digits */
char *z_to_str(Z a);                           /* arena string */
int z_cmp(Z a, Z b);
int z_cmpabs(Z a, Z b);
int z_is_zero(Z a);
int z_is_one(Z a);
int z_fits_i64(Z a, int64_t *out);
Z z_neg(Z a);
Z z_abs(Z a);
Z z_add(Z a, Z b);
Z z_sub(Z a, Z b);
Z z_mul(Z a, Z b);
Z z_mul_small(Z a, uint32_t m);
void z_divmod(Z a, Z b, Z *q, Z *r);           /* truncating: q toward zero, r has the sign of a */
Z z_div_round(Z a, Z b);                       /* nearest integer to a/b (ties away from zero) */
Z z_pow(Z a, unsigned e);
Z z_pow10(int64_t k);
Z z_mul_pow10(Z a, int64_t k);
int64_t z_digits(Z a);                         /* number of decimal digits of |a|, 0 for zero */
Z z_gcd(Z a, Z b);
Z z_iroot(Z a, unsigned r);                    /* floor of the r-th root of a >= 0 */
Z z_persist(Z a);                              /* copy into permanent memory */

/* ---- rationals ---- */
typedef struct { Z num, den; } Q;              /* den > 0, gcd(num, den) = 1 */

Q q_from_z(Z a);
Q q_make(Z num, Z den);
Q q_add(Q a, Q b);
Q q_sub(Q a, Q b);
Q q_mul(Q a, Q b);
Q q_div(Q a, Q b);
Q q_neg(Q a);
Q q_pow(Q a, int64_t e);
int q_sign(Q a);
int q_cmp(Q a, Q b);
int q_cmp_one(Q a);
int q_is_int(Q a);
char *q_to_str(Q a);
Q q_persist(Q a);

/* ---- polynomials in one letter, rational coefficients ---- */
typedef struct {
    int deg;                /* -1 for the zero polynomial */
    Q *c;                   /* c[i] is the coefficient of var^i */
    const char *var;        /* NULL for a constant */
} Poly;

Poly p_const(Q a);
Poly p_var(const char *name);
Poly p_add(Poly a, Poly b);
Poly p_sub(Poly a, Poly b);
Poly p_mul(Poly a, Poly b);
Poly p_pow(Poly a, unsigned e);
Poly p_scale(Poly a, Q s);
char *p_to_str(Poly a);

/* ---- quantities in letters: sums of k * a^e1 * b^e2 ..., rational k and rational exponents ---- */
#define NM_MAXL 12
typedef struct { Q k; Q e[NM_MAXL]; } CT;
typedef struct { int nt; CT *t; } C;

int letter_index(const char *name, size_t len);
const char *letter_name(int i);
int letter_count(void);
C c_zero(void);
C c_const(Q k);
C c_letter(int idx);
int c_is_zero(C a);
int c_const_value(C a, Q *out);
int c_uses(C a, int idx);
int c_is_monomial(C a);
C c_add(C a, C b);
C c_sub(C a, C b);
C c_neg(C a);
C c_scale(C a, Q k);
C c_mul(C a, C b);
C c_div(C a, C b);                            /* by a single term only */
C c_pow_int(C a, int64_t e);
int c_pow_q(C a, Q alpha, C *out);            /* 0 if not exact */
int c_equal(C a, C b);
C c_persist(C a);
C c_coeff_of(C a, int idx, Q e);              /* the part with letter idx to the power e, that letter removed */
char *c_to_str(C a);
char *ct_str(CT t, const char *extra_name, Q extra_e, int first);
char *c_term_str(C a, const char *var, Q e, int first);

/* ---- series: c[0] + c[1] x + ... + c[n-1] x^(n-1) + O(x^n), coefficients in letters ---- */
typedef struct { int n; C *c; } Ser;

Ser s_const(C a, int n);
Ser s_var(int n);
Ser s_add(Ser a, Ser b);
Ser s_sub(Ser a, Ser b);
Ser s_mul(Ser a, Ser b);
Ser s_div(Ser a, Ser b);
Ser s_scale(Ser a, C k);
Ser s_pow_int(Ser a, int64_t e);
Ser s_pow_q(Ser a, Q alpha);                  /* a(0) a single term with an exact power */
Ser s_deriv(Ser a);
Ser s_integ(Ser a);
Ser s_compose(Ser f, Ser g);                  /* f(g(x)), g(0) = 0 */
Ser s_trunc(Ser a, int n);
Ser s_persist(Ser a);
char *s_to_str(Ser a, const char *var, int last);   /* terms up to degree `last`, then O(...) */
int q_root_exact(Q a, int64_t n, Q *out);     /* a^(1/n) if rational */

/* ---- roots: an equation and a verified bracket ---- */
typedef struct Root {
    int deg;
    Z *c;                   /* integer coefficients, c[i] of y^i (permanent) */
    const char *var;
    Z X;                    /* current approximation X / 10^D (permanent) */
    int64_t D;
    int64_t w;              /* verified: a root lies strictly between (X - w)/10^D and (X + w)/10^D; 0 if exact */
    int certified;          /* 1 once a sign change confirmed the bracket at this D */
} Root;

char *parallelogram(C F, int xi, int yi, C *start, int have_start, int64_t order);

Root *root_new(Poly p, Q guess, int64_t start_places);
void root_refine(Root *r, int64_t places);    /* continue Newton's resolution to at least `places` decimals */
char *root_equation_str(Root *r);

/* ---- balls: midpoint m * 10^e, radius r * 10^e ---- */
typedef struct { Z m; int64_t e; Z r; } Ball;

Ball b_from_q(Q a, int64_t prec);
Ball b_from_root(Root *r, int64_t prec);
Ball b_add(Ball a, Ball b, int64_t prec);
Ball b_sub(Ball a, Ball b, int64_t prec);
Ball b_mul(Ball a, Ball b, int64_t prec);
Ball b_div(Ball a, Ball b, int64_t prec);
Ball b_pow(Ball a, int64_t e, int64_t prec);
int64_t b_guaranteed_places(Ball a, int64_t want);  /* largest k <= want with radius <= 0.5 * 10^-k, or -1 */
char *fixed_str(Z scaled, int64_t places);          /* scaled / 10^places written in decimal */
Z b_round_to_places(Ball a, int64_t places);

#endif
