/* Algebraic identities, primitive content, refusals, and an independent Q route. */
#include "nm.h"
#include <setjmp.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

extern jmp_buf nm_on_error;
extern char nm_error_msg[512];
static int checks, fail, identities, lowest, numeric, refusals;
static Q qi(int64_t n) { return q_from_z(z_from_i64(n)); }
static C N(int n) { return c_const(qi(n)); }
static C L(const char *s) { return c_letter(letter_index(s, strlen(s))); }
static void check(int ok, const char *label) {
    checks++;
    if (!ok) { fail++; printf("FAIL: %s\n", label); }
}
static void text_is(R r, const char *s) {
    lowest++; check(c_equal(poly_gcd(r.num,r.den),N(1)), "operation result in lowest terms");
    char *got = r_to_str(r);
    check(!strcmp(got, s), s);
    if (strcmp(got, s)) printf("  got: %s\n", got);
}
static void quotient(R r, C p, C q) {
    identities++; check(c_equal(c_mul(r.num, q), c_mul(p, r.den)), "quotient multiplied back");
    lowest++; check(c_equal(poly_gcd(r.num, r.den), N(1)), "lowest terms");
    check(q_sign(r.den.t[c_leading_index(r.den)].k) > 0, "positive leading denominator");
}
static Q eval_c(C p, Q *v) {
    Q s = qi(0);
    for (int t = 0; t < p.nt; t++) {
        Q m = p.t[t].k;
        for (int l = 0; l < letter_count(); l++) if (q_sign(p.t[t].e[l])) {
            int64_t e = 0;
            if (!z_fits_i64(p.t[t].e[l].num, &e)) nm_fail("test exponent too large");
            m = q_mul(m, q_pow(v[l], e));
        }
        s = q_add(s, m);
    }
    return s;
}
static uint32_t seed = 104729;
static Q sample(void) {
    seed = seed * 1664525u + 1013904223u;
    int n = (int)(seed % 19) - 9;
    seed = seed * 1664525u + 1013904223u;
    return q_make(z_from_i64(n), z_from_i64(seed % 5 + 1));
}
static void refuse(int which) {
    refusals++;
    jmp_buf saved; memcpy(saved, nm_on_error, sizeof saved);
    if (!setjmp(nm_on_error)) {
        C a = L("a"), x = L("x");
        if (which == 0) { C f; c_pow_q(x, q_make(z_from_i64(1), z_from_i64(2)), &f); (void)r_make(f, c_add(a, N(1))); }
        if (which == 1) (void)r_make(N(1), c_add(a, c_radical_q(qi(2), 2)));
        if (which == 2) (void)r_make(N(1), N(0));
        if (which == 3) (void)poly_exact_div(c_add(x, N(1)), x);
        check(0, "expected refusal");
    } else {
        const char *reasons[] = {"fractional exponents", "mixed with letters", "division by zero", "vitiose"};
        check(strstr(nm_error_msg, reasons[which]) != NULL, "refusal has the correct reason");
    }
    memcpy(nm_on_error, saved, sizeof saved);
}

int main(int argc, char **argv) {
    (void)argv;
    clock_t start = clock();
    if (setjmp(nm_on_error)) { printf("error: %s\n", nm_error_msg); return 1; }
    if (argc > 1) (void)L("x"); else (void)L("a");
    C a = L("a"), x = L("x"), b = L("b"), c = L("c"), d = L("d");
    int ia = letter_index("a",1), ix = letter_index("x",1), ib = letter_index("b",1);
    int ic = letter_index("c",1), id = letter_index("d",1);
    C apx = c_add(a,x), xma = c_sub(x,a), xx = c_pow_int(x,2), aa = c_pow_int(a,2);
    C p[11], q[11]; R r[11];
    p[0] = c_sub(xx,aa); q[0] = xma;
    p[1] = c_sub(c_pow_int(x,3),c_pow_int(a,3)); q[1] = xma;
    p[2] = c_sub(c_pow_int(x,4),c_pow_int(a,4)); q[2] = p[0];
    p[3] = c_sub(xx,N(1)); q[3] = c_pow_int(c_add(x,N(1)),2);
    p[4] = c_pow_int(apx,2); q[4] = p[0];
    p[5] = c_scale(x,qi(2)); q[5] = p[0];
    p[6] = N(1); q[6] = N(1);
    p[7] = c_mul(a,d); q[7] = c_mul(b,c);
    p[8] = c_mul(a,c_add(x,N(1))); q[8] = c_mul(b,c_add(x,N(1)));
    p[9] = aa; q[9] = c_add(b,x);
    p[10] = xma; q[10] = c_sub(aa,xx);
    for (int i = 0; i < 11; i++) r[i] = r_make(p[i],q[i]);
    r[5] = r_add(r_make(N(1),apx), r_make(N(1),xma));
    r[6] = r_add(r_make(a,c_sub(a,b)), r_make(b,c_sub(b,a)));
    r[7] = r_div(r_make(a,b), r_make(c,d));
    for (int i = 0; i < 11; i++) quotient(r[i],p[i],q[i]);
    const char *wanted[] = {"a + x", "a^2 + ax + x^2", "a^2 + x^2", "(x - 1)/(x + 1)",
        "(-a - x)/(a - x)", "-2x/(a^2 - x^2)", "1", "ad/bc", "a/b", "a^2/(b + x)", "-1/(a + x)"};
    for (int i = 0; i < 11; i++) text_is(r[i], wanted[i]);
    text_is(r_make(N(-1),apx), wanted[10]);
    text_is(r_make(apx,apx), "1"); text_is(r_make(N(0),apx), "0");
    text_is(r_make(N(2), c_scale(apx, qi(-2))), "-1/(a + x)");
    text_is(r_make(c_pow_int(x,-1),c_pow_int(a,-1)), "a/x");
    text_is(r_sub(r[5],r[5]), "0");
    text_is(r_mul(r[5],r_make(q[5],N(1))), "2x");
    text_is(r_pow_int(r_make(apx,xma),-2), "(a^2 - 2ax + x^2)/(a^2 + 2ax + x^2)");
    C f = c_mul(c_add(a,b),c_add(x,c)), g = c_mul(f,c_add(x,N(1)));
    C h = c_mul(f,c_add(x,N(2)));
    C gcd = poly_gcd(g,h);
    check(c_equal(gcd,f), "recursive nonconstant content in three letters");
    check(c_equal(poly_exact_div(g,gcd),c_add(x,N(1))), "gcd divides first input");
    check(c_equal(poly_exact_div(h,gcd),c_add(x,N(2))), "gcd divides second input");
    check(c_equal(poly_gcd(N(0),f),f), "gcd(0,p)");
    check(c_equal(poly_gcd(f,N(7)),N(1)), "gcd(p,constant)");
    check(c_is_zero(poly_gcd(N(0),N(0))), "gcd(0,0)");
    check(c_equal(poly_gcd(c_scale(c_mul(a,c_add(x,N(1))),qi(6)),c_scale(c_add(x,N(1)),qi(4))),c_add(x,N(1))), "answer-key gcd");
    check(c_equal(poly_gcd(c_pow_int(apx,3),c_pow_int(apx,2)),c_pow_int(apx,2)), "repeated common factor");
    for (int i = 0; i < 4; i++) refuse(i);
    R s = r_make(N(1), c_add(N(1),c_radical_q(qi(2),2)));
    check(c_equal(c_mul(s.num,c_add(N(1),c_radical_q(qi(2),2))),s.den), "existing surd inverse");
    for (int i = 0; i < 11; i++) r[i] = r_persist(r[i]);
    /* Original expressions below use Q alone, independently of polynomial/gcd code. */
    int sets = 0, skipped = 0;
    while (sets < 64) {
        arena_reset();
        Q v[NM_MAXL]; for (int l = 0; l < NM_MAXL; l++) v[l] = qi(0);
        Q A = sample(), X = sample(), B = sample(), Cc = sample(), D = sample();
        if (!q_sign(B) || !q_sign(Cc) || !q_sign(D) || !q_sign(q_sub(X,A)) || !q_sign(q_add(X,A)) ||
            !q_sign(q_add(X,qi(1))) || !q_sign(q_sub(A,B)) || !q_sign(q_add(B,X))) { skipped++; continue; }
        v[ia]=A; v[ix]=X; v[ib]=B; v[ic]=Cc; v[id]=D;
        Q want[11];
        want[0]=q_div(q_sub(q_pow(X,2),q_pow(A,2)),q_sub(X,A));
        want[1]=q_div(q_sub(q_pow(X,3),q_pow(A,3)),q_sub(X,A));
        want[2]=q_div(q_sub(q_pow(X,4),q_pow(A,4)),q_sub(q_pow(X,2),q_pow(A,2)));
        want[3]=q_div(q_sub(q_pow(X,2),qi(1)),q_add(q_add(q_pow(X,2),q_mul(qi(2),X)),qi(1)));
        want[4]=q_div(q_add(q_add(q_pow(X,2),q_mul(qi(2),q_mul(A,X))),q_pow(A,2)),q_sub(q_pow(X,2),q_pow(A,2)));
        want[5]=q_add(q_div(qi(1),q_add(X,A)),q_div(qi(1),q_sub(X,A)));
        want[6]=q_add(q_div(A,q_sub(A,B)),q_div(B,q_sub(B,A)));
        want[7]=q_div(q_div(A,B),q_div(Cc,D));
        want[8]=q_div(q_add(q_mul(A,X),A),q_add(q_mul(B,X),B));
        want[9]=q_div(q_pow(A,2),q_add(B,X));
        want[10]=q_div(q_sub(X,A),q_sub(q_pow(A,2),q_pow(X,2)));
        for (int i=0;i<11;i++) {
            Q got = q_div(eval_c(r[i].num,v),eval_c(r[i].den,v));
            numeric++; check(q_cmp(got,want[i])==0,"independent rational substitution");
        }
        sets++;
    }
    printf("ratfun: %d put-backs, %d lowest-terms, %d Q comparisons (%d sets, %d skipped), %d refusals; %.2f ms\n",
        identities,lowest,numeric,sets,skipped,refusals,1000.0*(clock()-start)/CLOCKS_PER_SEC);
    printf("unit ratfun: %d checks, %d failed\n",checks,fail);
    return fail != 0;
}
