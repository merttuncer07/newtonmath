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

struct Root;

/* ---- quantities in letters: sums of k * a^e1 * b^e2 ..., rational k and rational exponents ---- */
#define NM_MAXL 4096                          /* letters, surds and named values in one session */
#define CT_MAXV 12                             /* different letters in one term */
typedef struct { Q k; int n; short l[CT_MAXV]; Q x[CT_MAXV]; } CT;   /* k * product of l[i]^x[i]; l rising, x != 0 */
Q ct_e(const CT *t, int l);                    /* the exponent of letter l in the term (0 if absent) */
void ct_set(CT *t, int l, Q e);
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
C c_div(C a, C b);
C c_inv(C a);                                 /* a single term, or a quantity in surds (Euclid on the equation) */
int c_has_plain(C a);
C c_radical_q(Q c, int64_t n);                /* c^(1/n): exact, or a number times a surd */
C c_radical_c(C A, int64_t n);                /* A^(1/n) for numbers and surds: a surd over surds */
C c_imag_unit(void);
C c_surd_from_root(const char *name, struct Root *r);
int letter_is_surd(int i);
int letter_is_imag(int i);
struct Root *surd_root(int i);
int surd_radical(int i, C *A, int *n, int *neg);
C c_pow_int(C a, int64_t e);
int c_pow_q(C a, Q alpha, C *out);            /* 0 if not exact */
int c_equal(C a, C b);
C c_persist(C a);
C c_coeff_of(C a, int idx, Q e);              /* the part with letter idx to the power e, that letter removed */
char *c_to_str(C a);
int c_leading_index(C a);                    /* fixed degree/alphabetical term order; -1 for zero */
char *ct_str(CT t, const char *extra_name, Q extra_e, int first);
char *c_term_str(C a, const char *var, Q e, int first);

/* Polynomials in one letter with C coefficients, shared with elimination. */
typedef struct { int deg; C *c; } UP;
UP up_alloc(int deg);
UP up_trim(UP p);
UP up_from(C a, int v);
UP up_sub(UP a, UP b);

/* Rational functions over Q in ordinary letters. A unit denominator also holds existing C quantities. */
typedef struct { C num, den; } R;
C poly_gcd(C a, C b);                       /* monic gcd, including gcd(0,0) = 0 */
C poly_exact_div(C a, C b);                 /* refuses a nonzero remainder; checks by multiplication */
R r_from_c(C a);
R r_make(C num, C den);
R r_add(R a, R b);
R r_sub(R a, R b);
R r_mul(R a, R b);
R r_div(R a, R b);
R r_pow_int(R a, int64_t e);
int r_equal(R a, R b);
int r_is_zero(R a);
R r_persist(R a);
char *r_to_str(R a);

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

/* ---- systems: extermination by resultants (elim.c) ---- */
#define NM_MAXSOL 256
typedef struct {
    int nv, n;              /* unknowns; solutions found */
    C **val;                /* val[s][j]: the value of the j-th unknown in solution s */
    int curve;              /* some solutions are not isolated (a curve or more) */
    int unresolved;         /* roots not separated (complex roots of factors above degree two) */
    int rejected;           /* candidates brought in by extermination and refused by substitution */
} Solutions;
C det_nodiv(C *m, int n);
C elim_resultant(C A, C B, int v);
Solutions elim_solve(C *eqs, int ne, int *vars, int nv);
int elim_roots(C p, int v, C *out, int max, int *unresolved);

/* Exact conic factorization and closed-interval pole checks (elim.c). */
int elim_conic_factors(C p, int v, C *out, int max);
int elim_has_root_closed(C p, int v, Q lo, Q hi);

/* Rational integrals: a rational part and logarithmic/circular areas. */
enum { AREA_LOG, AREA_ATAN, AREA_ASIN };
typedef struct { R coef; C poly; int kind; } AreaTerm;
typedef struct { int var, n; R rational; AreaTerm *term; int root; int euler; C ek, ec; } Integral;
/* root: 1 + the letter put for sqrt(q), or 0; euler: 1 + the letter t = (sqrt(q) - ec)/(x - ek), or 0 */
typedef struct { int n; R *part; C *base; int *pw; } Apart;   /* part[i] over base[i]^pw[i]; pw 0: the whole part */
C integ_diff_poly(C p, int v);
C integ_subst_poly(C p, int v, C x);
R integ_derivative(Integral a);
Integral integ_rational(R f, int v);
Integral integ_persist(Integral a);
char *integ_to_str(Integral a);
C integ_value(Integral a, C x, C (*conic)(int, C));
Apart integ_apart(R f, int v);
Apart apart_persist(Apart a);
R apart_sum(Apart a);
char *apart_to_str(Apart a);
int letter_shown(const char *disp);

/* ---- areas under roots of the flowing letter (sqrtint.c) ---- */
Integral sqrt_integral(R f, int x, int s, C q);    /* f in x and s, where s = sqrt(q), q of degree 1 or 2 in x */
C sqrt_integral_value(Integral a, C q, C x, C (*conic)(int, C));
C asin_constant(C u, C (*conic)(int, C));
void sqrt_integral_split(R f, int s, C q, R *A, R *B);

/* ---- matrices (linalg.c) ---- */
typedef struct { int r, c; R *a; } Mat;       /* row by row */
Mat mat_new(int r, int c);
Mat mat_identity(int n);
Mat mat_add(Mat a, Mat b, int sign);
Mat mat_scale(Mat a, R k);
Mat mat_mul(Mat a, Mat b);
Mat mat_transpose(Mat a);
Mat mat_pow(Mat a, int64_t e);
Mat mat_inverse(Mat a);
Mat mat_nullspace(Mat a);
Mat mat_solve(Mat a, Mat b, Mat *nullspace);
int mat_equal(Mat a, Mat b);
int mat_rank(Mat a);
R mat_det(Mat a);
R *mat_charpoly(Mat a);                       /* det(t I - a): coefficients of t^0 .. t^n */
char *mat_to_str(Mat a);
char *mat_verdict(Mat input, Mat result, int rank); /* rank >= 0 includes the rank-drop locus */

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

/* ---- whole numbers (arith.c) ---- */
#define NM_MAXFAC 128
typedef struct { int sign, n; Z p[NM_MAXFAC]; int e[NM_MAXFAC]; int status[NM_MAXFAC]; } Factors;   /* status: 1 proved, 2 probable */
Z z_mod(Z a, Z m);
Z z_powmod(Z b, Z e, Z m);
Z z_xgcd(Z a, Z b, Z *s, Z *t);
int z_invmod(Z a, Z m, Z *out);
int z_crt(const Z *r, const Z *m, int k, Z *x, Z *M);
int z_isprime(Z n);                           /* 0 composite, 1 proved prime, 2 probable prime */
Factors z_factor(Z n);
Z z_sigma(Z n);
Z z_phi(Z n);
int z_divisors(Z n, Z *out, int max);
Z z_nextprime(Z n);

/* ---- complex balls and values from term rules (cplx.c) ---- */
typedef struct { Ball re, im; } CBall;
typedef struct { int s, T; Poly D; Poly *N; Poly h; } TermRule;   /* D(n) c_n = sum N_t(n) c_{n-t} - h_{n-s} */
CBall cb_from_q(Q re, Q im, int64_t prec);
CBall cb_add(CBall a, CBall b, int64_t prec);
CBall cb_sub(CBall a, CBall b, int64_t prec);
CBall cb_mul(CBall a, CBall b, int64_t prec);
CBall cb_div(CBall a, CBall b, int64_t prec);
Q cb_abs_upper(CBall a);
int rule_converges(const TermRule *R, Q r);
CBall rule_value(const TermRule *R, const Q *seed, int nseed, int64_t start, CBall A, int deriv, int64_t places, int64_t prec);
int letter_is_named(int i);
Ball named_ball(int i, int64_t prec);
C c_named(const char *disp, Ball (*value)(void *, int64_t), void *data);

#endif
