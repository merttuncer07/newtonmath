/* Areas under roots: each result is differentiated here again, by the test's own rule (s' = q'/(2s)), and
 * must give the curve back; arcs are checked through the square of their moment. Refusals give a reason. */
#include "nm.h"
#include <stdio.h>
#include <string.h>
#include <setjmp.h>

extern jmp_buf nm_on_error;
extern char nm_error_msg[512];
static int checks, failed, back;
static Q qi(int n) { return q_from_z(z_from_i64(n)); }
static C N(int n) { return c_const(qi(n)); }
static void ck(int ok, const char *s) { checks++; if (!ok) { failed++; printf("FAIL %s\n", s); } }
static int X, S;
static C q;

/* d/dx of a polynomial in x and s, as num/(2s) */
static C dnum(C p) { return c_add(c_mul(c_mul(N(2), c_letter(S)), integ_diff_poly(p, X)), c_mul(integ_diff_poly(q, X), integ_diff_poly(p, S))); }
static R d(R f) {
    C n = c_sub(c_mul(dnum(f.num), f.den), c_mul(f.num, dnum(f.den)));
    return r_make(n, c_mul(c_mul(N(2), c_letter(S)), c_mul(f.den, f.den)));
}
/* compare after s^2 -> q, by substituting s = q's root symbolically: multiply out with (a + b s) form */
static int same(R a, R b) {
    R A1, B1, A2, B2;
    sqrt_integral_split(a, S, q, &A1, &B1); sqrt_integral_split(b, S, q, &A2, &B2);
    return r_equal(A1, A2) && r_equal(B1, B2);
}
static void area(R f, const char *name) {
    Integral I = sqrt_integral(f, X, S, q);
    R g = d(I.rational);
    int arcs = 0;
    for (int i = 0; i < I.n; i++) {
        AreaTerm t = I.term[i];
        if (t.kind == AREA_LOG) g = r_add(g, r_mul(t.coef, d(r_make(t.poly, N(1))).num.nt ? r_div(d(r_make(t.poly, N(1))), r_make(t.poly, N(1))) : r_from_c(N(0))));
        else if (t.kind == AREA_ATAN) {
            C up = integ_diff_poly(t.poly, X);
            g = r_add(g, r_make(c_mul(c_div(t.coef.num, t.coef.den), up), c_add(N(1), c_mul(t.poly, t.poly))));
        } else arcs++;
    }
    if (!arcs) { back++; ck(same(g, f), name); return; }
    /* an arc c*asin(u): (f - g)^2 must be c^2 u'^2/(1 - u^2) */
    R rest = r_sub(f, g), want = r_from_c(N(0));
    for (int i = 0; i < I.n; i++) if (I.term[i].kind == AREA_ASIN) {
        C u = I.term[i].poly, up = integ_diff_poly(u, X);
        want = r_add(want, r_div(r_mul(r_mul(I.term[i].coef, I.term[i].coef), r_from_c(c_mul(up, up))), r_from_c(c_sub(N(1), c_mul(u, u)))));
    }
    back++; ck(same(r_mul(rest, rest), want), name);
}
static R subst(C p, int v, R T) {                   /* p with the letter v replaced by T */
    R r = r_from_c(c_zero());
    for (int i = 0; i < p.nt; i++) {
        C m; m.nt = 1; m.t = arena_alloc(sizeof(CT)); m.t[0] = p.t[i];
        int64_t e; z_fits_i64(p.t[i].e[v].num, &e); m.t[0].e[v] = qi(0);
        r = r_add(r, r_mul(r_from_c(m), r_pow_int(T, e)));
    }
    return r;
}
static void euler_back(R f) {
    Integral I = sqrt_integral(f, X, S, q);
    ck(I.euler != 0, "Euler's substitution used");
    int t = I.euler - 1;
    R T = r_div(r_sub(r_from_c(c_letter(S)), r_from_c(I.ec)), r_from_c(c_sub(c_letter(X), I.ek)));
    R g = d(I.rational);
    for (int i = 0; i < I.n; i++) {
        R P = subst(I.term[i].poly, t, T), dP = d(P);
        if (I.term[i].kind == AREA_LOG) g = r_add(g, r_mul(I.term[i].coef, r_div(dP, P)));
        else g = r_add(g, r_mul(I.term[i].coef, r_div(dP, r_add(r_from_c(N(1)), r_mul(P, P)))));
    }
    back++; ck(same(g, f), "Euler area put back");
}
static void refused(R f, const char *why) {
    jmp_buf saved; memcpy(saved, nm_on_error, sizeof saved);
    if (!setjmp(nm_on_error)) { (void)sqrt_integral(f, X, S, q); ck(0, why); }
    else ck(strstr(nm_error_msg, why) != NULL, why);
    memcpy(nm_on_error, saved, sizeof saved);
}
static void with(C rad) { q = rad; S = letter_shown(c_to_str(rad)); }

int main(void) {
    if (setjmp(nm_on_error)) { printf("error: %s (after %d checks)\n", nm_error_msg, checks); return 1; }
    X = letter_index("x", 1);
    C x = c_letter(X), x2 = c_mul(x, x);
    R s, one = r_from_c(N(1));
    /* hyperbolas */
    C hyps[] = {c_add(x2, N(1)), c_sub(x2, N(1)), c_add(c_add(x2, c_mul(N(2), x)), N(5)), c_add(c_mul(N(4), x2), N(1)), c_sub(x2, x)};
    for (int h = 0; h < 5; h++) {
        with(hyps[h]); s = r_from_c(c_letter(S));
        area(s, "sqrt(q)"); area(r_div(one, s), "1/sqrt(q)"); area(r_mul(r_from_c(x), s), "x sqrt(q)");
        area(r_mul(r_from_c(x2), s), "x^2 sqrt(q)"); area(r_div(r_from_c(c_add(c_mul(x2, x), N(3))), s), "(x^3 + 3)/sqrt(q)");
        area(r_add(r_div(one, r_from_c(c_add(x2, N(2)))), s), "rational beside the root");
    }
    /* circles */
    C circ[] = {c_sub(N(1), x2), c_sub(c_mul(N(2), x), x2), c_sub(N(3), c_mul(N(4), x2))};
    for (int h = 0; h < 3; h++) {
        with(circ[h]); s = r_from_c(c_letter(S));
        area(s, "circle sqrt(q)"); area(r_div(one, s), "circle 1/sqrt(q)"); area(r_div(r_from_c(c_mul(x2, x)), s), "circle x^3/sqrt(q)");
        area(r_mul(r_from_c(c_add(x, N(1))), s), "circle (x + 1) sqrt(q)");
    }
    /* straight lines: the root as the new letter */
    C lines[] = {c_add(x, N(1)), c_add(c_mul(N(2), x), N(1)), c_sub(N(4), c_mul(N(3), x))};   /* logs with rational roots: the test differentiates them as quotients */
    for (int h = 0; h < 3; h++) {
        with(lines[h]); s = r_from_c(c_letter(S));
        area(s, "line sqrt(q)"); area(r_div(one, r_mul(r_from_c(x), s)), "line 1/(x sqrt(q))");
        area(r_mul(r_from_c(x2), s), "line x^2 sqrt(q)"); area(r_div(one, r_add(one, s)), "line 1/(1 + sqrt(q))");
    }
    /* refusals */
    with(c_add(c_mul(x2, x), N(1))); refused(r_from_c(c_letter(S)), "degree 3");
    /* Euler's substitution: the result is brought back from t = (s - c)/(x - k) and differentiated */
    C eq[] = {c_add(x2, N(1)), c_add(x2, N(1)), c_sub(N(1), x2), c_add(x2, N(3))};
    C ed[] = {x, x2, x, c_add(x, N(1))};
    for (int h = 0; h < 4; h++) { with(eq[h]); euler_back(r_div(one, r_mul(r_from_c(ed[h]), r_from_c(c_letter(S))))); }
    with(c_sub(N(3), x2)); refused(r_div(one, r_mul(r_from_c(x), r_from_c(c_letter(S)))), "no rational point");
    with(c_add(c_add(x2, c_mul(N(2), x)), N(1))); refused(r_from_c(c_letter(S)), "perfect square");
    with(c_sub(N(-1), x2)); refused(r_from_c(c_letter(S)), "negative for every");
    with(c_add(x2, c_letter(letter_index("a", 1)))); refused(r_from_c(c_letter(S)), "letters under the root");
    printf("sqrtint: %d areas put back by the test's own moment\n", back);
    printf("unit sqrtint: %d checks, %d failed\n", checks, failed);
    return failed != 0;
}
