/* Module test for arith.c: residues, primes (proved or probable), factors, divisors, and the amicable numbers
 * of Newton's notebook (Descartes' rule). */
#include "nm.h"

#include <setjmp.h>
#include <stdio.h>
#include <string.h>

extern jmp_buf nm_on_error;
extern char nm_error_msg[512];

static Z D(const char *s) { return z_from_dec(s, strlen(s)); }
static Z zi(int64_t v) { return z_from_i64(v); }
static int fail, checks;
static void expect(const char *what, const char *got, const char *want) {
    checks++;
    int ok = !strcmp(got, want);
    printf("%-46s %s%s\n", what, got, ok ? "" : "   FAIL");
    if (!ok) { printf("   want %s\n", want); fail++; }
}
static char *fac(Factors F) {
    char *s = arena_alloc(4096), *o = s;
    *o = 0;
    for (int i = 0; i < F.n; i++) {
        o += sprintf(o, "%s%s", i ? " * " : "", z_to_str(F.p[i]));
        if (F.e[i] > 1) o += sprintf(o, "^%d", F.e[i]);
        if (F.status[i] == 2) o += sprintf(o, "(probable)");
    }
    return s;
}
static const char *verdict(int v) { return v == 1 ? "proved prime" : v == 2 ? "probable prime" : "composite"; }

int main(void) {
    if (setjmp(nm_on_error)) { printf("error: %s\n", nm_error_msg); return 1; }
    Z s, t;
    expect("gcd(240, 46)", z_to_str(z_xgcd(zi(240), zi(46), &s, &t)), "2");
    expect("   240 s + 46 t", z_to_str(z_add(z_mul(zi(240), s), z_mul(zi(46), t))), "2");
    Z inv; z_invmod(zi(17), zi(3120), &inv);
    expect("17^-1 mod 3120", z_to_str(inv), "2753");
    expect("2^1000 mod 1000007 (as bc gives)", z_to_str(z_powmod(zi(2), zi(1000), zi(1000007))), "783922");
    Z r[3] = {zi(2), zi(3), zi(2)}, m[3] = {zi(3), zi(5), zi(7)}, x, M;
    z_crt(r, m, 3, &x, &M);
    expect("x = 2 mod 3, 3 mod 5, 2 mod 7", z_to_str(x), "23");
    expect("561 (Carmichael)", verdict(z_isprime(zi(561))), "composite");
    expect("2^61 - 1", verdict(z_isprime(z_sub(z_pow(zi(2), 61), zi(1)))), "proved prime");
    expect("2^89 - 1 (beyond 3.3e24: Pocklington)", verdict(z_isprime(z_sub(z_pow(zi(2), 89), zi(1)))), "proved prime");
    expect("2^127 - 1", verdict(z_isprime(z_sub(z_pow(zi(2), 127), zi(1)))), "proved prime");
    expect("2^32 + 1 (Euler)", fac(z_factor(z_add(z_pow(zi(2), 32), zi(1)))), "641 * 6700417");
    expect("2^64 + 1 (Landry)", fac(z_factor(z_add(z_pow(zi(2), 64), zi(1)))), "274177 * 67280421310721");
    expect("10^20 + 1", fac(z_factor(z_add(z_pow(zi(10), 20), zi(1)))), "73 * 137 * 1676321 * 5964848081");
    expect("2^4 3^2 5 7^3 (whole factors)", fac(z_factor(zi(16 * 9 * 5 * 343))), "2^4 * 3^2 * 5 * 7^3");
    expect("phi(36)", z_to_str(z_phi(zi(36))), "12");
    Z dv[64]; int nd = z_divisors(zi(28), dv, 64);
    char buf[256] = ""; for (int i = 0; i < nd; i++) sprintf(buf + strlen(buf), "%s%s", i ? " " : "", z_to_str(dv[i]));
    expect("divisors of 28 (a perfect number)", buf, "1 2 4 7 14 28");
    expect("next prime after 10^12", z_to_str(z_nextprime(D("1000000000000"))), "1000000000039");
    /* Descartes' rule, as in Newton's notebook: if p = 3 2^(n-1) - 1, q = 3 2^n - 1, r = 9 2^(2n-1) - 1 are
     * prime, then 2^n p q and 2^n r are amicable: each is the sum of the other's aliquot parts */
    int ns[3] = {2, 4, 7};
    const char *want[3][2] = {{"220", "284"}, {"17296", "18416"}, {"9363584", "9437056"}};
    for (int k = 0; k < 3; k++) {
        int n = ns[k];
        Z p = z_sub(z_mul(zi(3), z_pow(zi(2), (unsigned)(n - 1))), zi(1));
        Z q = z_sub(z_mul(zi(3), z_pow(zi(2), (unsigned)n)), zi(1));
        Z rr = z_sub(z_mul(zi(9), z_pow(zi(2), (unsigned)(2 * n - 1))), zi(1));
        int primes = z_isprime(p) == 1 && z_isprime(q) == 1 && z_isprime(rr) == 1;
        Z A = z_mul(z_pow(zi(2), (unsigned)n), z_mul(p, q)), B = z_mul(z_pow(zi(2), (unsigned)n), rr);
        int amicable = !z_cmp(z_sub(z_sigma(A), A), B) && !z_cmp(z_sub(z_sigma(B), B), A);
        char line[128]; sprintf(line, "%s %s %s", z_to_str(A), z_to_str(B), primes && amicable ? "amicable" : "not");
        char wl[128]; sprintf(wl, "%s %s amicable", want[k][0], want[k][1]);
        char what[64]; sprintf(what, "Descartes' rule, n = %d", n);
        expect(what, line, wl);
    }
    printf("unit arith: %d checks, %d failed\n", checks, fail);
    return fail != 0;
}
