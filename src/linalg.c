/* Matrices of exact quantities: numbers, surds and letters.
 *
 * Linear equations are solved by exterminating one unknown after another, the way Newton taught them in his
 * algebra lectures (now called Gaussian elimination). Entries are rational functions; numbers and surds
 * keep unit denominators. Every inverse, solution and null vector is multiplied
 * back before it is returned. */
#include "nm.h"

#include <stdio.h>
#include <string.h>

static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }
static R rn(int64_t v) { return r_from_c(c_const(qi(v))); }
static R rneg(R a) { a.num = c_neg(a.num); return a; }
static R rinv(R a) { return r_make(a.den, a.num); }

Mat mat_new(int r, int c) {
    Mat m; m.r = r; m.c = c;
    m.a = arena_alloc((size_t)(r * c > 0 ? r * c : 1) * sizeof(R));
    for (int i = 0; i < r * c; i++) m.a[i] = rn(0);
    return m;
}

Mat mat_identity(int n) { Mat m = mat_new(n, n); for (int i = 0; i < n; i++) m.a[i * n + i] = rn(1); return m; }

#define AT(m, i, j) ((m).a[(i) * (m).c + (j)])

Mat mat_add(Mat a, Mat b, int sign) {
    if (a.r != b.r || a.c != b.c) nm_fail("matrices of sizes %dx%d and %dx%d cannot be added", a.r, a.c, b.r, b.c);
    Mat m = mat_new(a.r, a.c);
    for (int i = 0; i < a.r * a.c; i++) m.a[i] = sign > 0 ? r_add(a.a[i], b.a[i]) : r_sub(a.a[i], b.a[i]);
    return m;
}

Mat mat_scale(Mat a, R k) {
    Mat m = mat_new(a.r, a.c);
    for (int i = 0; i < a.r * a.c; i++) m.a[i] = r_mul(a.a[i], k);
    return m;
}

Mat mat_mul(Mat a, Mat b) {
    if (a.c != b.r) nm_fail("a %dx%d matrix times a %dx%d matrix: the inner sizes differ", a.r, a.c, b.r, b.c);
    Mat m = mat_new(a.r, b.c);
    for (int i = 0; i < a.r; i++)
        for (int j = 0; j < b.c; j++) {
            R s = rn(0);
            for (int k = 0; k < a.c; k++) if (!r_is_zero(AT(a, i, k)) && !r_is_zero(AT(b, k, j))) s = r_add(s, r_mul(AT(a, i, k), AT(b, k, j)));
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
    uint64_t n = e < 0 ? (uint64_t)(-(e + 1)) + 1 : (uint64_t)e;
    if (e < 0) a = mat_inverse(a);
    Mat r = mat_identity(a.r);
    while (n) { if (n & 1) r = mat_mul(r, a); n >>= 1; if (n) a = mat_mul(a, a); }
    return r;
}

int mat_equal(Mat a, Mat b) {
    if (a.r != b.r || a.c != b.c) return 0;
    for (int i = 0; i < a.r * a.c; i++) if (!r_equal(a.a[i], b.a[i])) return 0;
    return 1;
}

static int has_letters(Mat a) {
    for (int i = 0; i < a.r * a.c; i++) if (c_has_plain(a.a[i].num) || c_has_plain(a.a[i].den)) return 1;
    return 0;
}

/* ---------------- elimination over the field of rational functions ---------------- */

typedef struct { Mat m; int rank; int *pivcol; R det; } Echelon;

/* reduced row echelon form of a, with the determinant of its square part along the way */
static Echelon echelon(Mat a) {
    Echelon E;
    E.m = mat_new(a.r, a.c);
    memcpy(E.m.a, a.a, (size_t)(a.r * a.c) * sizeof(R));
    E.pivcol = arena_alloc((size_t)(a.r + 1) * sizeof(int));
    E.det = rn(1);
    int row = 0;
    for (int col = 0; col < a.c && row < a.r; col++) {
        int p = -1;
        for (int i = row; i < a.r; i++) if (!r_is_zero(AT(E.m, i, col))) { p = i; break; }
        if (p < 0) { E.det = rn(0); continue; }
        if (p != row) {                                /* exchange rows: the determinant changes sign */
            for (int j = 0; j < a.c; j++) { R t = AT(E.m, p, j); AT(E.m, p, j) = AT(E.m, row, j); AT(E.m, row, j) = t; }
            E.det = rneg(E.det);
        }
        R pv = AT(E.m, row, col);
        E.det = r_mul(E.det, pv);
        R inv = rinv(pv);
        for (int j = 0; j < a.c; j++) AT(E.m, row, j) = r_mul(AT(E.m, row, j), inv);
        for (int i = 0; i < a.r; i++) {                /* clear the column above and below */
            if (i == row || r_is_zero(AT(E.m, i, col))) continue;
            R f = AT(E.m, i, col);
            for (int j = 0; j < a.c; j++) AT(E.m, i, j) = r_sub(AT(E.m, i, j), r_mul(f, AT(E.m, row, j)));
        }
        E.pivcol[row++] = col;
    }
    E.rank = row;
    if (row < a.r) E.det = rn(0);
    return E;
}

R mat_det(Mat a) {
    if (a.r != a.c) nm_fail("only a square matrix has a determinant");
    if (has_letters(a)) {
        R *cp = mat_charpoly(a);
        return a.r % 2 ? rneg(cp[0]) : cp[0];
    }
    return echelon(a).det;
}

int mat_rank(Mat a) {
    return echelon(a).rank;
}

Mat mat_inverse(Mat a) {
    if (a.r != a.c) nm_fail("only a square matrix has an inverse");
    int n = a.r;
    Mat aug = mat_new(n, 2 * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) AT(aug, i, j) = AT(a, i, j);
        AT(aug, i, n + i) = rn(1);
    }
    Echelon E = echelon(aug);
    if (n && (E.rank < n || E.pivcol[n - 1] != n - 1)) nm_fail("the matrix is singular: its determinant is 0");
    Mat inv = mat_new(n, n);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) AT(inv, i, j) = AT(E.m, i, n + j);
    if (!mat_equal(mat_mul(a, inv), mat_identity(n))) nm_fail("internal check failed: A times its inverse is not 1 (vitiose)");
    return inv;
}

/* the null space: a basis of the vectors v with a v = 0, each checked */
Mat mat_nullspace(Mat a) {
    Echelon E = echelon(a);
    int free = a.c - E.rank;
    Mat N = mat_new(a.c, free);
    int *isp = arena_alloc((size_t)(a.c + 1) * sizeof(int));
    for (int j = 0; j < a.c; j++) isp[j] = -1;
    for (int r = 0; r < E.rank; r++) isp[E.pivcol[r]] = r;
    int k = 0;
    for (int f = 0; f < a.c; f++) {
        if (isp[f] >= 0) continue;
        AT(N, f, k) = rn(1);
        for (int r = 0; r < E.rank; r++) AT(N, E.pivcol[r], k) = rneg(AT(E.m, r, f));
        k++;
    }
    Mat z = mat_mul(a, N);
    for (int i = 0; i < z.r * z.c; i++) if (!r_is_zero(z.a[i])) nm_fail("internal check failed: a null vector is not null (vitiose)");
    return N;
}

/* a x = b: one solution x (a column), and the null space for the rest; fails if there is none */
Mat mat_solve(Mat a, Mat b, Mat *nullspace) {
    if (b.r != a.r) nm_fail("the right side has %d rows, the matrix %d", b.r, a.r);
    Mat aug = mat_new(a.r, a.c + b.c);
    for (int i = 0; i < a.r; i++) {
        for (int j = 0; j < a.c; j++) AT(aug, i, j) = AT(a, i, j);
        for (int j = 0; j < b.c; j++) AT(aug, i, a.c + j) = AT(b, i, j);
    }
    Echelon E = echelon(aug);
    for (int r = 0; r < E.rank; r++)
        if (E.pivcol[r] >= a.c) {
            if (has_letters(a) || has_letters(b)) nm_fail("for general values of the letters, the equations contradict each other; special values are not classified");
            nm_fail("the equations contradict each other: there is no solution");
        }
    Mat x = mat_new(a.c, b.c);
    for (int r = 0; r < E.rank; r++)
        for (int j = 0; j < b.c; j++) AT(x, E.pivcol[r], j) = AT(E.m, r, a.c + j);
    if (!mat_equal(mat_mul(a, x), b)) nm_fail("internal check failed: the solution does not satisfy the equations (vitiose)");
    if (nullspace) *nullspace = mat_nullspace(a);
    return x;
}

/* det(t I - a) as coefficients c[0..n] of t^0..t^n, by Berkowitz: sums and products only */
R *mat_charpoly(Mat a) {
    if (a.r != a.c) nm_fail("only a square matrix has a characteristic polynomial");
    int n = a.r;
    R *V = arena_alloc((size_t)(n + 2) * sizeof(R));
    V[0] = rn(1);
    int len = 1;
    if (n > 0) { V[1] = rneg(AT(a, 0, 0)); len = 2; }
    for (int r = 1; r < n; r++) {
        R *q = arena_alloc((size_t)(r + 2) * sizeof(R));
        q[0] = rn(1);
        q[1] = rneg(AT(a, r, r));
        R *w = arena_alloc((size_t)r * sizeof(R)), *w2 = arena_alloc((size_t)r * sizeof(R));
        for (int i = 0; i < r; i++) w[i] = AT(a, i, r);
        for (int k = 0; k < r; k++) {
            R s = rn(0);
            for (int i = 0; i < r; i++) s = r_add(s, r_mul(AT(a, r, i), w[i]));
            q[k + 2] = rneg(s);
            for (int i = 0; i < r; i++) {
                R t = rn(0);
                for (int j = 0; j < r; j++) t = r_add(t, r_mul(AT(a, i, j), w[j]));
                w2[i] = t;
            }
            R *tmp = w; w = w2; w2 = tmp;
        }
        R *NV = arena_alloc((size_t)(r + 3) * sizeof(R));
        for (int i = 0; i <= r + 1; i++) {
            R s = rn(0);
            for (int j = 0; j <= i && j < len; j++) s = r_add(s, r_mul(q[i - j], V[j]));
            NV[i] = s;
        }
        V = NV; len = r + 2;
    }
    R *c = arena_alloc((size_t)(n + 1) * sizeof(R));      /* V[k] is the coefficient of t^(n-k) */
    for (int k = 0; k <= n; k++) c[n - k] = V[k];
    return c;
}

char *mat_to_str(Mat a) {
    size_t cap = 16;
    char **cells = arena_alloc((size_t)(a.r * a.c + 1) * sizeof(char *));
    for (int i = 0; i < a.r * a.c; i++) { cells[i] = r_to_str(a.a[i]); cap += strlen(cells[i]) + 4; }
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

/* Rank drops exactly where every maximal nonzero minor vanishes, on the domain
 * of the input entries. A failed pivot alone is not a rank-drop condition. */
typedef struct Condition { C p; char *text; struct Condition *next; } Condition;
static void condition_add(Condition **head, C p) {
    if (c_is_zero(p)) return;
    p = c_scale(p, q_div(qi(1), p.t[c_leading_index(p)].k));
    for (Condition *c = *head; c; c = c->next) if (c_equal(c->p, p)) return;
    Condition *item = arena_alloc(sizeof *item);
    item->p = p; item->text = c_to_str(p); item->next = NULL;
    while (*head) head = &(*head)->next;
    *head = item;
}
static Condition *denominators(Mat m) {
    Condition *head = NULL;
    for (int i = 0; i < m.r * m.c; i++) if (c_has_plain(m.a[i].den)) condition_add(&head, m.a[i].den);
    return head;
}
static int next_combination(int *a, int k, int n) {
    for (int i = k - 1; i >= 0; i--) if (a[i] < n - k + i) {
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return 1;
    }
    return 0;
}
static Condition *rank_conditions(Mat m, int rank) {
    if (!rank) return NULL;
    int *rows = arena_alloc((size_t)rank * sizeof(int)), *cols = arena_alloc((size_t)rank * sizeof(int));
    for (int i = 0; i < rank; i++) rows[i] = i;
    Condition *head = NULL;
    do {
        for (int i = 0; i < rank; i++) cols[i] = i;
        do {
            Mat minor = mat_new(rank, rank);
            for (int i = 0; i < rank; i++) for (int j = 0; j < rank; j++) AT(minor,i,j) = AT(m,rows[i],cols[j]);
            R det = mat_det(minor);
            if (r_is_zero(det)) continue;
            if (!c_has_plain(det.num)) return NULL; /* a nonzero constant minor prevents a rank drop */
            condition_add(&head, det.num);
        } while (next_combination(cols, rank, m.c));
    } while (next_combination(rows, rank, m.r));
    return head;
}
static size_t condition_size(Condition *p) {
    size_t n = 0;
    for (; p; p = p->next) n += strlen(p->text) + 12;
    return n;
}
static char *conditions_print(char *out, Condition *p, const char *join) {
    int first = 1;
    for (; p; p = p->next) {
        out += sprintf(out, "%s%s = 0", first ? "" : join, p->text);
        first = 0;
    }
    return out;
}
char *mat_verdict(Mat input, Mat result, int rank) {
    if (!has_letters(input) && !has_letters(result)) {
        char *s = arena_alloc(8); strcpy(s,"[exact]"); return s;
    }
    int order[NM_MAXL], n = 0;
    size_t cap = 384;
    for (int l = 0; l < letter_count(); l++) {
        if (letter_is_surd(l) || letter_is_named(l)) continue;
        int used = 0;
        Mat ms[2] = {input,result};
        for (int k = 0; k < 2; k++) for (int i = 0; i < ms[k].r * ms[k].c; i++)
            if (c_uses(ms[k].a[i].num,l) || c_uses(ms[k].a[i].den,l)) used = 1;
        if (used) {
            int i = n++;
            while (i > 0 && strcmp(letter_name(l),letter_name(order[i-1])) < 0) { order[i] = order[i-1]; i--; }
            order[i] = l; cap += strlen(letter_name(l)) + 2;
        }
    }
    Condition *domain = denominators(input), *poles = denominators(result);
    Condition *drop = rank >= 0 ? rank_conditions(input,rank) : NULL;
    cap += condition_size(domain) + condition_size(poles) + condition_size(drop);
    char *s = arena_alloc(cap), *o = s;
    o += sprintf(o,"[exact, for general ");
    for (int i = 0; i < n; i++) o += sprintf(o,"%s%s",i ? ", " : "",letter_name(order[i]));
    if (domain) { o += sprintf(o,"; input undefined where "); o = conditions_print(o,domain," or "); }
    if (rank >= 0) {
        if (drop) { o += sprintf(o,"; rank drops where "); o = conditions_print(o,drop," and "); }
        else o += sprintf(o,"; rank never drops on the input domain");
    }
    if (poles) { o += sprintf(o,"; displayed formula undefined where "); o = conditions_print(o,poles," or "); }
    sprintf(o,"]");
    return s;
}
