/* Module test for cplx.c: values of functions from their term rules at complex points and of their fluxions,
 * against bc -l (the expected digits below were computed with bc, rounded to 30 places). */
#include "nm.h"

#include <setjmp.h>
#include <stdio.h>
#include <string.h>

extern jmp_buf nm_on_error;
extern char nm_error_msg[512];

static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }
static Q qf(int64_t a, int64_t b) { return q_make(z_from_i64(a), z_from_i64(b)); }
static Poly P1(int deg, const int64_t *c) {       /* a polynomial in n with whole coefficients */
    Poly p; p.deg = deg; p.var = "n"; p.c = arena_alloc((size_t)(deg + 1 > 0 ? deg + 1 : 1) * sizeof(Q));
    for (int i = 0; i <= deg; i++) p.c[i] = qi(c[i]);
    return p;
}
static int fail, checks;
static TermRule Eg; static Q egseed[1];
static Ball exp_at_1(void *data, int64_t prec) {      /* the value behind the letter exp(1) */
    (void)data;
    return rule_value(&Eg, egseed, 1, 1, cb_from_q(q_from_z(z_from_i64(1)), q_from_z(z_zero()), prec), 0, prec - 20, prec).re;
}

static void expect(const char *what, Ball b, const char *want) {
    checks++;
    char *got = fixed_str(b_round_to_places(b, 30), 30);
    char digits[128]; int k = 0;
    for (const char *p = got; *p; p++) if (*p != '.') digits[k++] = *p;
    digits[k] = 0;
    char *d = digits;
    int neg = d[0] == '-'; if (neg) d++;
    while (*d == '0' && d[1]) d++;
    char norm[130]; snprintf(norm, sizeof norm, "%s%s", neg ? "-" : "", d);
    int ok = !strcmp(norm, want) && b_guaranteed_places(b, 30) == 30;
    printf("%-40s %s%s\n", what, got, ok ? "" : "   FAIL");
    if (!ok) { printf("   want %s\n", want); fail++; }
}

int main(void) {
    if (setjmp(nm_on_error)) { printf("error: %s\n", nm_error_msg); return 1; }
    int64_t P = 60, W = 30;
    /* exp: n c_n = c_{n-1} */
    TermRule E; E.s = 1; E.T = 1; int64_t eD[2] = {0, 1}, eN[1] = {1};
    E.D = P1(1, eD); E.N = arena_alloc(2 * sizeof(Poly)); E.N[1] = P1(0, eN); E.h = P1(-1, NULL);
    Q eseed[1] = {qi(1)};
    /* sin: n(n-1) c_n = -c_{n-2} */
    TermRule S; S.s = 2; S.T = 2; int64_t sD[3] = {0, -1, 1}, sN2[1] = {-1};
    S.D = P1(2, sD); S.N = arena_alloc(3 * sizeof(Poly)); S.N[1] = P1(-1, NULL); S.N[2] = P1(0, sN2); S.h = P1(-1, NULL);
    Q sseed[2] = {qi(0), qi(1)};
    /* log(1 + x): n c_n = -(n - 1) c_{n-1} - h, h = -1 */
    TermRule G; G.s = 1; G.T = 1; int64_t gD[2] = {0, 1}, gN[2] = {1, -1}, gh[1] = {-1};
    G.D = P1(1, gD); G.N = arena_alloc(2 * sizeof(Poly)); G.N[1] = P1(1, gN); G.h = P1(0, gh);
    Q gseed[2] = {qi(0), qi(1)};

    CBall i1 = cb_from_q(qi(0), qi(1), P);
    CBall v = rule_value(&E, eseed, 1, 1, i1, 0, W, P);
    expect("exp(i): real part = cos 1", v.re, "540302305868139717400936607443");
    expect("        imaginary part = sin 1", v.im, "841470984807896506652502321630");
    v = rule_value(&E, eseed, 1, 1, cb_from_q(qi(1), qi(2), P), 0, W, P);
    expect("exp(1 + 2i): real part", v.re, "-1131204383756813638431255255511");
    expect("            imaginary part", v.im, "2471726672004818927616930893552");
    v = rule_value(&S, sseed, 2, 2, i1, 0, W, P);
    expect("sin(i) = i sinh 1: imaginary part", v.im, "1175201193643801456882381850596");
    expect("                   real part", v.re, "0");
    v = rule_value(&S, sseed, 2, 2, cb_from_q(qf(1, 2), qi(0), P), 1, W, P);
    expect("sin'(1/2) = cos(1/2)", v.re, "877582561890372716116281582604");
    v = rule_value(&E, eseed, 1, 1, cb_from_q(qi(1), qi(0), P), 2, W, P);
    expect("exp''(1) = e", v.re, "2718281828459045235360287471353");
    v = rule_value(&G, gseed, 2, 2, cb_from_q(qi(0), qf(1, 2), P), 0, W, P);
    expect("log(1 + i/2): real = log(sqrt(5)/2)", v.re, "111571775657104877883147545155");
    expect("              imaginary = atan(1/2)", v.im, "463647609000806116214256231461");
    /* on the edge: refused, with the reason */
    checks++;
    jmp_buf saved; memcpy(saved, nm_on_error, sizeof saved);
    if (!setjmp(nm_on_error)) { rule_value(&G, gseed, 2, 2, i1, 0, W, P); printf("log(1 + i) on the edge not refused   FAIL\n"); fail++; }
    else printf("log(1 + i): %s\n", nm_error_msg);
    memcpy(nm_on_error, saved, sizeof saved);
    /* a named value: exp(1) is a letter; its square stays exact, its value comes from the rule */
    Eg = E; egseed[0] = qi(1);
    C e1 = c_named("exp(1)", exp_at_1, NULL);
    C e2 = c_mul(e1, e1);
    checks++;
    printf("exp(1)^2 as a quantity: %s%s\n", c_to_str(e2), strcmp(c_to_str(e2), "exp(1)^2") ? "   FAIL" : "");
    if (strcmp(c_to_str(e2), "exp(1)^2")) fail++;
    int li = -1;
    for (int l = 0; l < letter_count(); l++) if (letter_is_named(l)) li = l;
    Ball b = named_ball(li, P);
    expect("the value of exp(1)", b, "2718281828459045235360287471353");
    printf("unit cplx: %d checks, %d failed\n", checks, fail);
    return fail != 0;
}
