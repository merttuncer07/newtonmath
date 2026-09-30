/* Matrices of exact quantities: numbers, surds and letters.
 *
 * Linear equations are solved by exterminating one unknown after another, the way Newton taught them in his
 * algebra lectures (now called Gaussian elimination). Where letters stand in the entries, a quotient of two
 * sums of letters would be needed, so determinant and characteristic polynomial are taken without any division
 * (Berkowitz) and the rest waits for rational functions. Every inverse, solution and null vector is multiplied
 * back before it is returned. */
#include "nm.h"

#include <stdio.h>
#include <string.h>

static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }

Mat mat_new(int r, int c) {
    Mat m; m.r = r; m.c = c;
    m.a = arena_alloc((size_t)(r * c > 0 ? r * c : 1) * sizeof(C));
    for (int i = 0; i < r * c; i++) m.a[i] = c_zero();
    return m;
}

Mat mat_identity(int n) { Mat m = mat_new(n, n); for (int i = 0; i < n; i++) m.a[i * n + i] = c_const(qi(1)); return m; }

#define AT(m, i, j) ((m).a[(i) * (m).c + (j)])

Mat mat_add(Mat a, Mat b, int sign) {
    if (a.r != b.r || a.c != b.c) nm_fail("matrices of sizes %dx%d and %dx%d cannot be added", a.r, a.c, b.r, b.c);
    Mat m = mat_new(a.r, a.c);
    for (int i = 0; i < a.r * a.c; i++) m.a[i] = sign > 0 ? c_add(a.a[i], b.a[i]) : c_sub(a.a[i], b.a[i]);
    return m;
}

Mat mat_scale(Mat a, C k) {
    Mat m = mat_new(a.r, a.c);
    for (int i = 0; i < a.r * a.c; i++) m.a[i] = c_mul(a.a[i], k);
    return m;
}

Mat mat_mul(Mat a, Mat b) {
    if (a.c != b.r) nm_fail("a %dx%d matrix times a %dx%d matrix: the inner sizes differ", a.r, a.c, b.r, b.c);
    Mat m = mat_new(a.r, b.c);
    for (int i = 0; i < a.r; i++)
        for (int j = 0; j < b.c; j++) {
            C s = c_zero();
            for (int k = 0; k < a.c; k++) if (!c_is_zero(AT(a, i, k)) && !c_is_zero(AT(b, k, j))) s = c_add(s, c_mul(AT(a, i, k), AT(b, k, j)));
            AT(m, i, j) = s;
        }
    return m;
}

Mat mat_transpose(Mat a) {
    Mat m = mat_new(a.c, a.r);
    for (int i = 0; i < a.r; i++) for (int j = 0; j < a.c; j++) AT(m, j, i) = AT(a, i, j);
    return m;
}

Mat mat_pow(Mat a, int64_t e) {
    if (a.r != a.c) nm_fail("only a square matrix has powers");
    if (e < 0) return mat_pow(mat_inverse(a), -e);
    Mat r = mat_identity(a.r);
    while (e) { if (e & 1) r = mat_mul(r, a); e >>= 1; if (e) a = mat_mul(a, a); }
    return r;
}

int mat_equal(Mat a, Mat b) {
    if (a.r != b.r || a.c != b.c) return 0;
    for (int i = 0; i < a.r * a.c; i++) if (!c_equal(a.a[i], b.a[i])) return 0;
    return 1;
}

static int has_letters(Mat a) {
    for (int i = 0; i < a.r * a.c; i++) if (c_has_plain(a.a[i])) return 1;
    return 0;
}

/* ---------------- elimination over numbers and surds ---------------- */

typedef struct { Mat m; int rank; int *pivcol; C det; } Echelon;

/* reduced row echelon form of a (numbers and surds), with the determinant of its square part along the way */
static Echelon echelon(Mat a) {
    Echelon E;
    E.m = mat_new(a.r, a.c);
    memcpy(E.m.a, a.a, (size_t)(a.r * a.c) * sizeof(C));
    E.pivcol = arena_alloc((size_t)(a.r + 1) * sizeof(int));
    E.det = c_const(qi(1));
    int row = 0;
    for (int col = 0; col < a.c && row < a.r; col++) {
        int p = -1;
        for (int i = row; i < a.r; i++) if (!c_is_zero(AT(E.m, i, col))) { p = i; break; }
        if (p < 0) { E.det = c_zero(); continue; }
        if (p != row) {                                /* exchange rows: the determinant changes sign */
            for (int j = 0; j < a.c; j++) { C t = AT(E.m, p, j); AT(E.m, p, j) = AT(E.m, row, j); AT(E.m, row, j) = t; }
            E.det = c_neg(E.det);
        }
        C pv = AT(E.m, row, col);
        E.det = c_mul(E.det, pv);
        C inv = c_inv(pv);
        for (int j = 0; j < a.c; j++) AT(E.m, row, j) = c_mul(AT(E.m, row, j), inv);
        for (int i = 0; i < a.r; i++) {                /* clear the column above and below */
            if (i == row || c_is_zero(AT(E.m, i, col))) continue;
            C f = AT(E.m, i, col);
            for (int j = 0; j < a.c; j++) AT(E.m, i, j) = c_sub(AT(E.m, i, j), c_mul(f, AT(E.m, row, j)));
        }
        E.pivcol[row++] = col;
    }
    E.rank = row;
    if (row < a.r) E.det = c_zero();
    return E;
}

C mat_det(Mat a) {
    if (a.r != a.c) nm_fail("only a square matrix has a determinant");
    if (has_letters(a)) return det_nodiv(a.a, a.r);
    return echelon(a).det;
}

int mat_rank(Mat a) {
    if (has_letters(a)) nm_fail("the rank of a matrix with letters depends on their values (this comes later)");
    return echelon(a).rank;
}

Mat mat_inverse(Mat a) {
    if (a.r != a.c) nm_fail("only a square matrix has an inverse");
    if (has_letters(a)) nm_fail("the inverse of a matrix with letters needs quotients of sums of letters (this comes later)");
    int n = a.r;
    Mat aug = mat_new(n, 2 * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) AT(aug, i, j) = AT(a, i, j);
        AT(aug, i, n + i) = c_const(qi(1));
    }
    Echelon E = echelon(aug);
    if (E.rank < n || E.pivcol[n - 1] != n - 1) nm_fail("the matrix is singular: its determinant is 0");
    Mat inv = mat_new(n, n);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) AT(inv, i, j) = AT(E.m, i, n + j);
    if (!mat_equal(mat_mul(a, inv), mat_identity(n))) nm_fail("internal check failed: A times its inverse is not 1 (vitiose)");
    return inv;
}

/* the null space: a basis of the vectors v with a v = 0, each checked */
Mat mat_nullspace(Mat a) {
    if (has_letters(a)) nm_fail("the null space of a matrix with letters depends on their values (this comes later)");
    Echelon E = echelon(a);
    int free = a.c - E.rank;
    Mat N = mat_new(a.c, free);
    int *isp = arena_alloc((size_t)(a.c + 1) * sizeof(int));
    for (int j = 0; j < a.c; j++) isp[j] = -1;
    for (int r = 0; r < E.rank; r++) isp[E.pivcol[r]] = r;
    int k = 0;
    for (int f = 0; f < a.c; f++) {
        if (isp[f] >= 0) continue;
        AT(N, f, k) = c_const(qi(1));
        for (int r = 0; r < E.rank; r++) AT(N, E.pivcol[r], k) = c_neg(AT(E.m, r, f));
        k++;
    }
    Mat z = mat_mul(a, N);
    for (int i = 0; i < z.r * z.c; i++) if (!c_is_zero(z.a[i])) nm_fail("internal check failed: a null vector is not null (vitiose)");
    return N;
}

/* a x = b: one solution x (a column), and the null space for the rest; fails if there is none */
Mat mat_solve(Mat a, Mat b, Mat *nullspace) {
    if (b.r != a.r) nm_fail("the right side has %d rows, the matrix %d", b.r, a.r);
    if (has_letters(a) || has_letters(b)) nm_fail("solving with letters needs quotients of sums of letters (this comes later)");
    Mat aug = mat_new(a.r, a.c + b.c);
    for (int i = 0; i < a.r; i++) {
        for (int j = 0; j < a.c; j++) AT(aug, i, j) = AT(a, i, j);
        for (int j = 0; j < b.c; j++) AT(aug, i, a.c + j) = AT(b, i, j);
    }
    Echelon E = echelon(aug);
    for (int r = 0; r < E.rank; r++)
        if (E.pivcol[r] >= a.c) nm_fail("the equations contradict each other: there is no solution");
    Mat x = mat_new(a.c, b.c);
    for (int r = 0; r < E.rank; r++)
        for (int j = 0; j < b.c; j++) AT(x, E.pivcol[r], j) = AT(E.m, r, a.c + j);
    if (!mat_equal(mat_mul(a, x), b)) nm_fail("internal check failed: the solution does not satisfy the equations (vitiose)");
    if (nullspace) *nullspace = mat_nullspace(a);
    return x;
}

/* det(t I - a) as coefficients c[0..n] of t^0..t^n, by Berkowitz: sums and products only */
C *mat_charpoly(Mat a) {
    if (a.r != a.c) nm_fail("only a square matrix has a characteristic polynomial");
    int n = a.r;
    C *V = arena_alloc((size_t)(n + 2) * sizeof(C));
    V[0] = c_const(qi(1));
    int len = 1;
    if (n > 0) { V[1] = c_neg(AT(a, 0, 0)); len = 2; }
    for (int r = 1; r < n; r++) {
        C *q = arena_alloc((size_t)(r + 2) * sizeof(C));
        q[0] = c_const(qi(1));
        q[1] = c_neg(AT(a, r, r));
        C *w = arena_alloc((size_t)r * sizeof(C)), *w2 = arena_alloc((size_t)r * sizeof(C));
        for (int i = 0; i < r; i++) w[i] = AT(a, i, r);
        for (int k = 0; k < r; k++) {
            C s = c_zero();
            for (int i = 0; i < r; i++) s = c_add(s, c_mul(AT(a, r, i), w[i]));
            q[k + 2] = c_neg(s);
            for (int i = 0; i < r; i++) {
                C t = c_zero();
                for (int j = 0; j < r; j++) t = c_add(t, c_mul(AT(a, i, j), w[j]));
                w2[i] = t;
            }
            C *tmp = w; w = w2; w2 = tmp;
        }
        C *NV = arena_alloc((size_t)(r + 3) * sizeof(C));
        for (int i = 0; i <= r + 1; i++) {
            C s = c_zero();
            for (int j = 0; j <= i && j < len; j++) s = c_add(s, c_mul(q[i - j], V[j]));
            NV[i] = s;
        }
        V = NV; len = r + 2;
    }
    C *c = arena_alloc((size_t)(n + 1) * sizeof(C));      /* V[k] is the coefficient of t^(n-k) */
    for (int k = 0; k <= n; k++) c[n - k] = V[k];
    return c;
}

char *mat_to_str(Mat a) {
    size_t cap = 16;
    char **cells = arena_alloc((size_t)(a.r * a.c + 1) * sizeof(char *));
    for (int i = 0; i < a.r * a.c; i++) { cells[i] = c_to_str(a.a[i]); cap += strlen(cells[i]) + 4; }
    char *s = arena_alloc(cap + (size_t)a.r * 4), *o = s;
    *o++ = '[';
    for (int i = 0; i < a.r; i++) {
        o += sprintf(o, "%s[", i ? ", " : "");
        for (int j = 0; j < a.c; j++) o += sprintf(o, "%s%s", j ? ", " : "", cells[i * a.c + j]);
        *o++ = ']';
    }
    *o++ = ']';
    *o = 0;
    return s;
}
