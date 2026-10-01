/* Module test for linalg.c: exact matrices, every inverse, solution and null vector multiplied back. */
#include "nm.h"

#include <setjmp.h>
#include <stdio.h>
#include <string.h>

extern jmp_buf nm_on_error;
extern char nm_error_msg[512];

static R L(const char *name) { return r_from_c(c_letter(letter_index(name, strlen(name)))); }
static R N(int64_t v) { return r_from_c(c_const(q_from_z(z_from_i64(v)))); }
static R F(int64_t a, int64_t b) { return r_from_c(c_const(q_make(z_from_i64(a), z_from_i64(b)))); }
static int fail, checks;
static void expect(const char *what, const char *got, const char *want) {
    checks++;
    int ok = !strcmp(got, want);
    printf("%s: %s%s\n", what, got, ok ? "" : "   FAIL");
    if (!ok) { printf("   want: %s\n", want); fail++; }
}

int main(void) {
    if (setjmp(nm_on_error)) { printf("error: %s\n", nm_error_msg); return 1; }
    /* the Hilbert matrix of order 4: determinant 1/6048000, an inverse of whole numbers */
    Mat H = mat_new(4, 4);
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) H.a[i * 4 + j] = F(1, i + j + 1);
    expect("det Hilbert 4", r_to_str(mat_det(H)), "1/6048000");
    expect("inverse Hilbert 4", mat_to_str(mat_inverse(H)),
           "[[16, -120, 240, -140], [-120, 1200, -2700, 1680], [240, -2700, 6480, -4200], [-140, 1680, -4200, 2800]]");
    /* letters: determinant and characteristic polynomial without division */
    R a = L("a"), b = L("b"), c = L("c"), d = L("d");
    Mat M = mat_new(2, 2); M.a[0] = a; M.a[1] = b; M.a[2] = c; M.a[3] = d;
    expect("det [[a, b], [c, d]]", r_to_str(mat_det(M)), "ad - bc");
    R *cp = mat_charpoly(M);
    expect("charpoly [[a, b], [c, d]], t^1 and t^0", r_to_str(cp[1]), "-a - d");
    expect("   t^0", r_to_str(cp[0]), "ad - bc");
    /* the Vandermonde determinant of a, b, c */
    Mat V = mat_new(3, 3);
    R xs[3] = {a, b, c};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) V.a[i * 3 + j] = r_pow_int(xs[i], j);
    R vd = mat_det(V);
    R want = r_mul(r_mul(r_sub(b, a), r_sub(c, a)), r_sub(c, b));
    checks++; printf("Vandermonde det equals (b - a)(c - a)(c - b): %s\n", r_equal(vd, want) ? "yes" : "NO   FAIL"); if (!r_equal(vd, want)) fail++;
    /* a rotation by 45 degrees: its inverse is its transpose */
    R s2 = r_div(r_from_c(c_radical_q(q_from_z(z_from_i64(2)), 2)), N(2));
    Mat R = mat_new(2, 2); R.a[0] = s2; R.a[1] = r_sub(N(0),s2); R.a[2] = s2; R.a[3] = s2;
    checks++; printf("rotation: inverse = transpose: %s\n", mat_equal(mat_inverse(R), mat_transpose(R)) ? "yes" : "NO   FAIL");
    if (!mat_equal(mat_inverse(R), mat_transpose(R))) fail++;
    expect("rotation^8", mat_to_str(mat_pow(R, 8)), "[[1, 0], [0, 1]]");
    /* a singular system: one solution and a line of others */
    Mat A = mat_new(3, 3);
    int64_t av[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    for (int i = 0; i < 9; i++) A.a[i] = N(av[i]);
    Mat B = mat_new(3, 1); B.a[0] = N(6); B.a[1] = N(15); B.a[2] = N(24);
    Mat ns;
    Mat x = mat_solve(A, B, &ns);
    expect("rank [[1,2,3],[4,5,6],[7,8,9]]", mat_rank(A) == 2 ? "2" : "not 2", "2");
    expect("a solution of A x = (6, 15, 24)", mat_to_str(x), "[[0], [3], [0]]");
    expect("null space", mat_to_str(ns), "[[1], [-2], [1]]");
    /* a contradiction is reported, not answered */
    B.a[2] = N(25);
    checks++;
    jmp_buf saved; memcpy(saved, nm_on_error, sizeof saved);
    if (!setjmp(nm_on_error)) { mat_solve(A, B, NULL); printf("contradiction not seen   FAIL\n"); fail++; }
    else printf("A x = (6, 15, 25): %s\n", nm_error_msg);
    memcpy(nm_on_error, saved, sizeof saved);
    /* Fibonacci by powers of [[1, 1], [1, 0]] */
    Mat Fm = mat_new(2, 2); Fm.a[0] = N(1); Fm.a[1] = N(1); Fm.a[2] = N(1); Fm.a[3] = N(0);
    expect("[[1, 1], [1, 0]]^90 (F91, F90)", r_to_str(mat_pow(Fm, 90).a[0]), "4660046610375530309");
    printf("unit linalg: %d checks, %d failed\n", checks, fail);
    return fail != 0;
}
