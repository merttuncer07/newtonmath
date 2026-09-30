/* Module test for elim.c: systems solved by extermination, every solution checked by substitution. */
#include "nm.h"

#include <setjmp.h>
#include <stdio.h>
#include <string.h>

extern jmp_buf nm_on_error;
extern char nm_error_msg[512];

static C L(const char *name) { return c_letter(letter_index(name, strlen(name))); }
static C N(int64_t v) { return c_const(q_from_z(z_from_i64(v))); }
static C add(C a, C b) { return c_add(a, b); }
static C sub(C a, C b) { return c_sub(a, b); }
static C mul(C a, C b) { return c_mul(a, b); }
static C pw(C a, int e) { return c_pow_int(a, e); }

static int fail, checks;

static void show(const char *title, Solutions S, int expect_n) {
    printf("%s: %d solution%s", title, S.n, S.n == 1 ? "" : "s");
    if (S.rejected) printf(" (%d refused by substitution)", S.rejected);
    if (S.unresolved) printf(" (%d roots not separated)", S.unresolved);
    if (S.curve) printf(" (not isolated)");
    printf("\n");
    for (int s = 0; s < S.n; s++) {
        printf("   ");
        for (int j = 0; j < S.nv; j++) printf("%s%s", j ? ",  " : "", c_to_str(S.val[s][j]));
        printf("\n");
    }
    checks++;
    if (expect_n >= 0 && S.n != expect_n) { printf("   FAIL: expected %d\n", expect_n); fail++; }
}

int main(void) {
    if (setjmp(nm_on_error)) { printf("error: %s\n", nm_error_msg); return 1; }
    C x = L("x"), y = L("y"), z = L("z"), a = L("a");
    int ix = letter_index("x", 1), iy = letter_index("y", 1), iz = letter_index("z", 1);
    int xy[2] = {ix, iy}, xyz[3] = {ix, iy, iz};

    /* a circle and a line */
    C e1[2] = {sub(add(pw(x, 2), pw(y, 2)), N(1)), sub(y, x)};
    show("x^2 + y^2 = 1, y = x", elim_solve(e1, 2, xy, 2), 2);

    /* Newton's two curves meeting: y^2 = x^3 and y = x */
    C e2[2] = {sub(pw(y, 2), pw(x, 3)), sub(y, x)};
    show("y^2 = x^3, y = x", elim_solve(e2, 2, xy, 2), 2);

    /* two conics: four points */
    C e3[2] = {sub(add(pw(x, 2), pw(y, 2)), N(5)), sub(mul(x, y), N(2))};
    show("x^2 + y^2 = 5, xy = 2", elim_solve(e3, 2, xy, 2), 4);

    /* the elementary symmetric functions of 1, 2, 3: six solutions */
    C e4[3] = {sub(add(add(x, y), z), N(6)), sub(add(add(mul(x, y), mul(y, z)), mul(z, x)), N(11)), sub(mul(mul(x, y), z), N(6))};
    show("x + y + z = 6, xy + yz + zx = 11, xyz = 6", elim_solve(e4, 3, xyz, 3), 6);

    /* a cubic that does not factor: its real root as a surd; two complex roots left unresolved */
    C e5[2] = {sub(sub(pw(x, 3), x), N(1)), sub(y, pw(x, 2))};
    show("x^3 - x - 1 = 0, y = x^2", elim_solve(e5, 2, xy, 2), 1);

    /* no real point: a circle and a line apart; the solutions are complex */
    C e6[2] = {sub(add(pw(x, 2), pw(y, 2)), N(1)), sub(y, N(2))};
    show("x^2 + y^2 = 1, y = 2", elim_solve(e6, 2, xy, 2), 2);

    /* a letter in the coefficients: the line y = a x meets the parabola y = x^2 */
    C e7[2] = {sub(y, mul(a, x)), sub(y, pw(x, 2))};
    Solutions s7 = elim_solve(e7, 2, xy, 2);
    show("y = ax, y = x^2  (the letter a stays)", s7, -1);

    /* the same line twice: a curve of solutions */
    C e8[2] = {sub(y, x), sub(mul(N(2), y), mul(N(2), x))};
    show("y = x, 2y = 2x", elim_solve(e8, 2, xy, 2), 0);

    /* the resultant alone: exterminating y from y^2 = x and y = x - 2 */
    C r = elim_resultant(sub(pw(y, 2), x), sub(y, sub(x, N(2))), iy);
    printf("resultant of y^2 - x and y - x + 2 in y: %s\n", c_to_str(r));
    checks++;
    if (!c_equal(r, add(sub(pw(x, 2), mul(N(5), x)), N(4)))) { printf("   FAIL\n"); fail++; }

    printf("unit elim: %d checks, %d failed\n", checks, fail);
    return fail != 0;
}
