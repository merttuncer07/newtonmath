/* The language: reading a statement, working it out, and writing the result.
 *
 *   statement := "use" NAME | "let" NAME "=" expr [to] | expr [to]
 *   to        := "to" INTEGER "places" | "to" LETTER "^" INTEGER
 *   expr      := term { ("+" | "-") term }
 *   term      := unary { ("*" | "/") unary | unary }          (side by side multiplies: 2y, 3(x+1))
 *   unary     := "-" unary | power
 *   power     := primary [ "^" unary ]
 *   primary   := NUMBER | NAME{'} | NAME "(" expr ")" | "(" expr ")" | "sqrt" "(" expr ")"
 *              | "d/d"LETTER primary | "integral" "(" expr ["," LETTER] ")"
 *              | "root" "of" expr "=" expr ( "near" expr | { "," NAME{'} "(" expr ")" "=" expr } )
 *
 * One engine does the hard work: resolution. A number is resolved by Newton's iteration on its places (approx.c);
 * a series is resolved term by term, the next term found from the lowest terms left over, as in the Methodus.
 * The same `root of` serves an algebraic equation and a fluxional one. Plain ASCII is canonical; √ × · − ² ³ are
 * read as aliases. */
#include "nm.h"

#include <ctype.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern jmp_buf nm_on_error;
extern char nm_error_msg[512];

#define DEFAULT_PLACES 20
#define DEFAULT_ORDER 8
#define GUARD_DIGITS 30

/* ---------------- tokens ---------------- */

enum { T_END, T_NUM, T_NAME, T_OP, T_KEY };
typedef struct { int kind; const char *s; size_t len; char op; int primes; const char *at; } Tok;   /* at: where in the line */

static const char *KEYWORDS[] = {"let", "to", "places", "root", "of", "near", "sqrt", "integral", "use", "for", "starting", "sum", "if", "otherwise", NULL};

static Tok *toks;
static int ntok, pos;

static int is_key(const char *s, size_t len) {
    for (int i = 0; KEYWORDS[i]; i++)
        if (strlen(KEYWORDS[i]) == len && strncmp(KEYWORDS[i], s, len) == 0) return 1;
    return 0;
}

static void push(int kind, const char *s, size_t len, char op) {
    Tok t = {kind, s, len, op, 0, s};
    toks[ntok++] = t;
}

static void lex(const char *src) {
    size_t n = strlen(src);
    toks = arena_alloc((n * 3 + 2) * sizeof(Tok));
    ntok = 0; pos = 0;
    for (size_t i = 0; i < n;) {
        unsigned char c = (unsigned char)src[i];
        if (c == '#') break;
        if (isspace(c)) { i++; continue; }
        if (isdigit(c) || (c == '.' && i + 1 < n && isdigit((unsigned char)src[i + 1]))) {
            size_t j = i;
            while (j < n && isdigit((unsigned char)src[j])) j++;
            if (j < n && src[j] == '.') { j++; while (j < n && isdigit((unsigned char)src[j])) j++; }
            push(T_NUM, src + i, j - i, 0);
            i = j; continue;
        }
        if (isalpha(c) || c == '_') {
            size_t j = i;
            while (j < n && (isalnum((unsigned char)src[j]) || src[j] == '_')) j++;
            push(is_key(src + i, j - i) ? T_KEY : T_NAME, src + i, j - i, 0);
            while (j < n && src[j] == '\'') { toks[ntok - 1].primes++; j++; }   /* y', y'' */
            i = j; continue;
        }
        if ((c == '<' || c == '>') && i + 1 < n && src[i + 1] == '=') {
            push(T_OP, src + i, 2, c == '<' ? 'l' : 'g'); i += 2; continue;
        }
        if (strchr("+-*/^()=,[]<>", c)) { push(T_OP, src + i, 1, (char)c); i++; continue; }
        if (!strncmp(src + i, "\xE2\x88\x9A", 3)) { push(T_KEY, "sqrt", 4, 0); toks[ntok - 1].at = src + i; i += 3; continue; }      /* √ */
        if (!strncmp(src + i, "\xC3\x97", 2) || !strncmp(src + i, "\xC2\xB7", 2)) {                     /* × · */
            push(T_OP, "*", 1, '*'); toks[ntok - 1].at = src + i; i += 2; continue;
        }
        if (!strncmp(src + i, "\xE2\x88\x92", 3)) { push(T_OP, "-", 1, '-'); toks[ntok - 1].at = src + i; i += 3; continue; }         /* − */
        if (!strncmp(src + i, "\xC2\xB2", 2) || !strncmp(src + i, "\xC2\xB3", 2)) {                     /* ² ³ */
            push(T_OP, "^", 1, '^'); toks[ntok - 1].at = src + i;
            push(T_NUM, src[i + 1] == '\xB2' ? "2" : "3", 1, 0); toks[ntok - 1].at = src + i + 1;
            i += 2; continue;
        }
        nm_fail("unexpected character '%c'", c);
    }
    push(T_END, src + n, 0, 0);
}

static Tok *peek(void) { return &toks[pos]; }
static Tok *peek_at(int k) { return pos + k < ntok ? &toks[pos + k] : &toks[ntok - 1]; }
static int at_op(char op) { return peek()->kind == T_OP && peek()->op == op; }
static int at_key(const char *k) {
    return peek()->kind == T_KEY && strlen(k) == peek()->len && !strncmp(peek()->s, k, peek()->len);
}
static void expect_op(char op) {
    if (!at_op(op)) nm_fail("expected '%c'", op);
    pos++;
}
static void expect_key(const char *k) {
    if (!at_key(k)) nm_fail("expected '%s'", k);
    pos++;
}

/* ---------------- syntax tree ---------------- */

enum { N_NUM, N_NAME, N_NEG, N_BIN, N_SQRT, N_ROOT, N_APPLY, N_DERIV, N_INTEG, N_INDEX, N_SUM, N_CASES, N_CMP };

typedef struct Cond { const char *s; size_t len; int primes; struct Node *at, *val; } Cond;

typedef struct Node {
    int kind; char op;
    const char *s; size_t len; int primes;      /* a name, or the letter of d/dx and integral */
    struct Node *a, *b, *c;
    Cond *conds; int nconds;
    int para; const char *us; size_t ul; struct Node *start;   /* root of ... for y [starting y = ...] */
    struct Node **args; int nargs;                               /* f(x, y) */
    struct Node **cval, **ccond; int ncases;                     /* v1 if c1, v2 if c2, ..., v otherwise */
} Node;

static int same_text(const char *a, size_t la, const char *b, size_t lb) { return la == lb && !strncmp(a, b, la); }

/* Newton's letters: "initialibus literis a, b, c" the given quantities, "finalibus literis v, x, y, z" the flowing */
static int flowing_letter(const char *s, size_t len) { return len == 1 && strchr("tuvwxz", s[0]); }

static Node *mk(int kind) { Node *n = arena_alloc(sizeof *n); memset(n, 0, sizeof *n); n->kind = kind; return n; }
static Node *bin(char op, Node *a, Node *b) { Node *n = mk(N_BIN); n->op = op; n->a = a; n->b = b; return n; }

static Node *expr(void);
static Node *unary(void);
static Node *power(void);

static Node *primary(void) {
    Tok *t = peek();
    if (t->kind == T_NUM) { Node *n = mk(N_NUM); n->s = t->s; n->len = t->len; pos++; return n; }
    if (t->kind == T_NAME) {
        Tok *t1 = peek_at(1), *t2 = peek_at(2);
        if (t->len == 1 && t->s[0] == 'd' && t1->kind == T_OP && t1->op == '/' && t2->kind == T_NAME
            && t2->len >= 2 && t2->s[0] == 'd') {                      /* d/dx */
            pos += 3;
            Node *n = mk(N_DERIV);
            n->s = t2->s + 1; n->len = t2->len - 1;
            if (at_op('(')) { pos++; n->a = expr(); expect_op(')'); }
            else n->a = power();
            return n;
        }
        pos++;
        if (at_op('[')) {                                              /* A[k]: a term of a sequence */
            pos++;
            Node *n = mk(N_INDEX);
            n->s = t->s; n->len = t->len;
            n->a = expr(); expect_op(']');
            return n;
        }
        if (at_op('(')) {                                              /* f(x), f(x, y), or y (x + 1) */
            pos++;
            Node *n = mk(N_APPLY);
            n->s = t->s; n->len = t->len; n->primes = t->primes;
            n->args = arena_alloc(16 * sizeof(Node *));
            n->args[n->nargs++] = n->a = expr();
            while (at_op(',')) {
                pos++;
                if (n->nargs == 16) nm_fail("too many arguments");
                n->args[n->nargs++] = expr();
            }
            expect_op(')');
            return n;
        }
        Node *n = mk(N_NAME); n->s = t->s; n->len = t->len; n->primes = t->primes;
        return n;
    }
    if (at_op('(')) { pos++; Node *e = expr(); expect_op(')'); return e; }
    if (at_key("sqrt")) {
        pos++;
        Node *n = mk(N_SQRT);
        if (at_op('(')) { pos++; n->a = expr(); expect_op(')'); }
        else n->a = primary();                          /* √2 */
        return n;
    }
    if (at_key("sum")) {                                               /* sum(A[k] x^k for k = 0 to 8) */
        pos++; expect_op('(');
        Node *n = mk(N_SUM);
        n->a = expr();
        expect_key("for");
        if (peek()->kind != T_NAME) nm_fail("expected the letter that runs, as in 'for k = 0 to 8'");
        n->s = peek()->s; n->len = peek()->len; pos++;
        expect_op('='); n->b = expr(); expect_key("to"); n->c = expr(); expect_op(')');
        return n;
    }
    if (at_key("integral")) {
        pos++; expect_op('(');
        Node *n = mk(N_INTEG);
        n->a = expr();
        if (at_op(',')) {
            pos++;
            if (peek()->kind != T_NAME) nm_fail("expected the letter to integrate in");
            n->s = peek()->s; n->len = peek()->len; pos++;
        }
        expect_op(')');
        return n;
    }
    if (at_key("root")) {
        pos++; expect_key("of");
        Node *n = mk(N_ROOT);
        n->a = expr(); expect_op('='); n->b = expr();
        if (at_key("near")) { pos++; n->c = expr(); return n; }
        if (at_key("for") || at_key("starting")) {                 /* Newton's parallelogram */
            n->para = 1;
            if (at_key("for")) {
                pos++;
                if (peek()->kind != T_NAME) nm_fail("expected the unknown after 'for'");
                n->us = peek()->s; n->ul = peek()->len; pos++;
            }
            if (at_key("starting")) {
                pos++;
                if (peek()->kind != T_NAME) nm_fail("expected 'starting y = ...'");
                if (n->us && !same_text(peek()->s, peek()->len, n->us, n->ul)) nm_fail("the start must be for %.*s", (int)n->ul, n->us);
                n->us = peek()->s; n->ul = peek()->len; pos++;
                expect_op('=');
                n->start = expr();
            }
            return n;
        }
        n->conds = arena_alloc(16 * sizeof(Cond));
        while (at_op(',')) {
            pos++;
            if (n->nconds == 16) nm_fail("too many starting conditions");
            Cond *c = &n->conds[n->nconds++];
            if (peek()->kind != T_NAME) nm_fail("expected a start such as y(0) = 1");
            c->s = peek()->s; c->len = peek()->len; c->primes = peek()->primes; pos++;
            expect_op('('); c->at = expr(); expect_op(')'); expect_op('='); c->val = expr();
        }
        if (!n->nconds) nm_fail("give a starting value: 'near 2' for a number, or 'y(0) = 1' for a series");
        return n;
    }
    if (t->kind == T_END) nm_fail("the statement ends too early");
    nm_fail("unexpected '%.*s'", (int)t->len, t->s);
    return NULL;
}

static Node *power(void) {
    Node *b = primary();
    if (at_op('^')) { pos++; return bin('^', b, unary()); }
    return b;
}

static Node *unary(void) {
    if (at_op('-')) { pos++; Node *n = mk(N_NEG); n->a = unary(); return n; }
    if (at_op('+')) { pos++; return unary(); }
    return power();
}

static Node *term(void) {
    Node *a = unary();
    for (;;) {
        if (at_op('*') || at_op('/')) { char op = peek()->op; pos++; a = bin(op, a, unary()); }
        else if (peek()->kind == T_NAME || peek()->kind == T_NUM || at_op('(') || at_key("sqrt") || at_key("root")
                 || at_key("integral"))
            a = bin('*', a, power());                   /* side by side: 2y, 3(x+1) */
        else return a;
    }
}

static Node *expr(void) {
    Node *a = term();
    while (at_op('+') || at_op('-')) { char op = peek()->op; pos++; a = bin(op, a, term()); }
    return a;
}

/* a condition: a = b, a < b, a > b, a <= b, a >= b */
static Node *cond(void) {
    Node *n = mk(N_CMP);
    n->a = expr();
    if (!(at_op('=') || at_op('<') || at_op('>') || at_op('l') || at_op('g'))) nm_fail("expected a comparison (=, <, >, <=, >=)");
    n->op = peek()->op; pos++;
    n->b = expr();
    return n;
}

/* a value, or values by cases: v1 if c1, v2 if c2, ..., v otherwise */
static Node *cases(void) {
    Node *e = expr();
    if (!at_key("if")) return e;
    Node *n = mk(N_CASES);
    n->cval = arena_alloc(32 * sizeof(Node *)); n->ccond = arena_alloc(32 * sizeof(Node *));
    pos++;
    n->cval[0] = e; n->ccond[0] = cond(); n->ncases = 1;
    while (at_op(',')) {
        pos++;
        Node *v = expr();
        if (at_key("otherwise")) { pos++; n->c = v; return n; }
        expect_key("if");
        if (n->ncases == 32) nm_fail("too many cases");
        n->cval[n->ncases] = v; n->ccond[n->ncases] = cond(); n->ncases++;
    }
    nm_fail("a definition by cases ends with '..., value otherwise'");
    return NULL;
}

/* parse a stored definition without disturbing the statement being read */
static Node *parse_text(const char *src) {
    Tok *st = toks; int sn = ntok, sp = pos;
    lex(src);
    Node *e = cases();
    if (peek()->kind != T_END) nm_fail("internal: stored definition did not parse");
    toks = st; ntok = sn; pos = sp;
    return e;
}

/* ---------------- values ---------------- */

enum { V_Q, V_POLY, V_ROOT, V_BALL, V_REC, V_NUMREC, V_FUNC, V_SEQ };   /* recipes: series, number; rules */
struct Def;
typedef struct { int kind; Q q; C c; Root *root; Ball ball; const char *src, *var; struct Def *def; } Val;   /* V_POLY holds c */

/* A quantity with its moment: a + b·o, where o·o is rejected (Methodus, Problem 1). The moment carries the
 * derivative with respect to the unknown term, which is what resolution needs. */
typedef struct { Ser a, b; int hasb; } Dual;

/* a rule f(x, y) = body, or a sequence A[0] = ..., A[k] = body, whose terms are kept in a table built forward */
typedef struct Def {
    int np; char **pn;                          /* parameters, or the index letter of a sequence */
    char *src;                                  /* the body */
    int ninit; int64_t *iidx; char **isrc;      /* starting terms of a sequence */
    int gen; int known; int cap; Val *table;    /* the table, valid while no name has been redefined */
    long building; int64_t progress;
} Def;

typedef struct Binding { char *name; Val v; struct Binding *next; } Binding;
static Binding *names;

static Binding *lookup(const char *s, size_t len) {
    for (Binding *b = names; b; b = b->next)
        if (strlen(b->name) == len && !strncmp(b->name, s, len)) return b;
    return NULL;
}

static int same(const char *a, size_t la, const char *b, size_t lb) { return la == lb && !strncmp(a, b, la); }

static char *dup_perm(const char *s, size_t len) { char *r = perm_alloc(len + 1); memcpy(r, s, len); r[len] = 0; return r; }

static Val vq(Q q) { Val v; memset(&v, 0, sizeof v); v.kind = V_Q; v.q = q; return v; }
static Val vc(C c) {
    Q k;
    if (c_const_value(c, &k)) return vq(k);
    Val v; memset(&v, 0, sizeof v); v.kind = V_POLY; v.c = c; return v;
}
static Val vball(Ball b) { Val v; memset(&v, 0, sizeof v); v.kind = V_BALL; v.ball = b; return v; }

static int64_t work_prec;                               /* significant digits for balls in this statement */

typedef struct { const char *s; size_t len; int hasv; Val v; int hasd; Dual d; } Param;
static Param *pstack;
static int pdepth, pcap, frame_base, calldepth;
static int let_gen;                                     /* bumped by every definition: tables start again */
static long stmt_id;

static Param *param_find(const char *s, size_t len) {
    for (int i = pdepth - 1; i >= frame_base; i--)
        if (pstack[i].len == len && !strncmp(pstack[i].s, s, len)) return &pstack[i];
    return NULL;
}

static Param *param_push(const char *s, size_t len) {
    if (pdepth == pcap) {
        int ncap = pcap ? pcap * 2 : 256;
        Param *np = perm_alloc((size_t)ncap * sizeof(Param));
        if (pdepth) memcpy(np, pstack, (size_t)pdepth * sizeof(Param));
        pstack = np; pcap = ncap;                      /* the old stack is left behind: small, and rare */
    }
    Param *p = &pstack[pdepth++];
    memset(p, 0, sizeof *p);
    p->s = s; p->len = len;
    return p;
}
static jmp_buf need_series;                             /* the exact pass meets something only a series can hold */

static C as_c(Val v) { return v.kind == V_POLY ? v.c : c_const(v.q); }

static Ball ball_of_c(C a, int64_t prec);

static Ball as_ball(Val v) {
    switch (v.kind) {
    case V_POLY: return ball_of_c(v.c, work_prec);
    case V_Q: return b_from_q(v.q, work_prec);
    case V_ROOT: return b_from_root(v.root, work_prec);
    case V_BALL: return v.ball;
    default: nm_fail("an unknown letter cannot be mixed with approximate numbers");
    }
    return v.ball;
}

static int is_exact(Val v) { return v.kind == V_Q || v.kind == V_POLY; }

/* the n-th root of every number between lo > 0 and hi, as one ball: roots of both ends by Newton's resolution */
static Ball ball_root_between(Q lo, Q hi, int64_t n, int64_t prec) {
    Root *ends[2];
    Q qs[2] = {lo, hi};
    for (int k = 0; k < 2; k++) {
        Poly p; p.deg = (int)n; p.var = "y";
        p.c = arena_alloc((size_t)(n + 1) * sizeof(Q));
        for (int i = 0; i <= n; i++) p.c[i] = q_from_z(z_zero());
        p.c[n] = q_from_z(z_from_i64(1)); p.c[0] = q_neg(qs[k]);
        Z g = z_iroot(z_div_round(z_mul_pow10(qs[k].num, 6 * n), qs[k].den), (unsigned)n);
        ends[k] = root_new(p, q_make(g, z_pow10(6)), 6);
        root_refine(ends[k], prec + 5);
        if (!ends[k]->certified) nm_fail("internal: a root of a radicand could not be certified");
    }
    int64_t D = ends[0]->D > ends[1]->D ? ends[0]->D : ends[1]->D;
    Z L = z_sub(z_mul_pow10(ends[0]->X, D - ends[0]->D), z_mul_pow10(z_from_i64(ends[0]->w), D - ends[0]->D));
    Z U = z_add(z_mul_pow10(ends[1]->X, D - ends[1]->D), z_mul_pow10(z_from_i64(ends[1]->w), D - ends[1]->D));
    Ball b;
    b.m = z_div_round(z_add(L, U), z_from_i64(2));
    b.r = z_add(z_div_round(z_sub(U, L), z_from_i64(2)), z_from_i64(1));
    b.e = -D;
    return b;
}

static Ball ball_of_surd(int l, int64_t prec) {
    if (letter_is_imag(l)) nm_fail("a complex number has no single real value here (functions at complex numbers come later)");
    Root *r = surd_root(l);
    if (r) return b_from_root(r, prec);
    C A; int n, neg;
    if (!surd_radical(l, &A, &n, &neg)) nm_fail("internal: a surd with no way to its value");
    Ball a = ball_of_c(A, prec + 10);
    Z lo = z_sub(a.m, a.r), hi = z_add(a.m, a.r);
    if (lo.s <= 0) nm_fail("the sign of %s could not be fixed at this precision", c_to_str(A));
    Q den = a.e >= 0 ? q_from_z(z_pow10(a.e)) : q_make(z_from_i64(1), z_pow10(-a.e));
    return ball_root_between(q_mul(q_from_z(lo), den), q_mul(q_from_z(hi), den), n, prec);
}

static Ball ball_of_c(C a, int64_t prec) {
    if (c_has_plain(a)) nm_fail("an unknown letter cannot be mixed with approximate numbers");
    Ball acc = b_from_q(q_from_z(z_zero()), prec);
    for (int t = 0; t < a.nt; t++) {
        Ball term = b_from_q(a.t[t].k, prec);
        for (int l = 0; l < letter_count(); l++) {
            int64_t e;
            if (!q_sign(a.t[t].e[l])) continue;
            if (!z_fits_i64(a.t[t].e[l].num, &e) || !q_is_int(a.t[t].e[l])) nm_fail("internal: a fractional power of a surd");
            term = b_mul(term, b_pow(ball_of_surd(l, prec), e, prec), prec);
        }
        acc = b_add(acc, term, prec);
    }
    return acc;
}

/* a quantity with i, as its real and imaginary parts */
static int split_imag(C a, C *re, C *im) {
    int li = -1;
    for (int l = 0; l < letter_count(); l++) if (letter_is_imag(l) && c_uses(a, l)) li = l;
    if (li < 0) return 0;
    *re = c_coeff_of(a, li, q_from_z(z_zero()));
    *im = c_coeff_of(a, li, q_from_z(z_from_i64(1)));
    return 1;
}
static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }

static Q parse_number(const char *s, size_t len) {
    const char *dot = memchr(s, '.', len);
    if (!dot) return q_from_z(z_from_dec(s, len));
    size_t ip = (size_t)(dot - s), fp = len - ip - 1;
    char *all = arena_alloc(len + 2);
    memcpy(all, s, ip); memcpy(all + ip, dot + 1, fp);
    if (ip + fp == 0) { all[0] = '0'; ip = 1; }
    return q_make(z_from_dec(all, ip + fp), z_pow10((int64_t)fp));   /* 0.1 is exactly 1/10 */
}

/* c^(1/n) for a rational c: exact, or a number times a surd (a letter with its equation) */
static Val vc(C c);
static Val nth_root(Q c, int64_t n) {
    if (n < 2 || n > 1000) nm_fail("root index out of range");
    return vc(c_radical_q(c, n));
}

/* a quantity in one letter with whole powers, as a polynomial (var: the letter required, or NULL for any) */
static Poly c_to_poly(C c, const char *var) {
    int li = -1;
    for (int l = 0; l < letter_count(); l++)
        if (c_uses(c, l)) {
            if (li >= 0 || (var && strcmp(letter_name(l), var))) {
                if (var) nm_fail("the equation's coefficients must be polynomials in %s alone", var);
                nm_fail("an equation for a number must have one unknown letter");
            }
            li = l;
        }
    int deg = 0;
    for (int t = 0; t < c.nt; t++) {
        int64_t e = 0;
        if (li >= 0 && (!q_is_int(c.t[t].e[li]) || !z_fits_i64(c.t[t].e[li].num, &e) || e < 0 || e > 100000))
            nm_fail("the equation must have whole powers of %s", letter_name(li));
        if (e > deg) deg = (int)e;
    }
    Poly p;
    p.var = li >= 0 ? letter_name(li) : (var ? var : NULL);
    p.deg = c.nt ? deg : -1;
    p.c = arena_alloc((size_t)(deg + 1) * sizeof(Q));
    for (int i = 0; i <= deg; i++) p.c[i] = qi(0);
    for (int t = 0; t < c.nt; t++) {
        int64_t e = 0;
        if (li >= 0) z_fits_i64(c.t[t].e[li].num, &e);
        p.c[e] = q_add(p.c[e], c.t[t].k);
    }
    while (p.deg >= 0 && q_sign(p.c[p.deg]) == 0) p.deg--;
    return p;
}

/* ---------------- the exact pass: numbers, rationals, polynomials ---------------- */

static Val eval(Node *n);
static Val series_value(Binding *b, Val arg);

/* while reading an equation's coefficients, the unknown and its derivatives stand for given numbers */
static struct { int active; const char *u; size_t ul; Q vals[8]; } ovr;

static Val power_val(Val base, Val ex) {
    if (ex.kind != V_Q) nm_fail("the exponent must be an exact number");
    Q e = ex.q;
    int64_t num, den;
    if (!z_fits_i64(e.num, &num) || !z_fits_i64(e.den, &den) || num > 1000000 || num < -1000000)
        nm_fail("exponent too large");
    if (base.kind == V_Q) {
        if (den == 1) return vq(q_pow(base.q, num));
        return nth_root(q_pow(base.q, num), den);
    }
    if (base.kind == V_POLY) {
        if (den == 1 && (num >= 0 || c_is_monomial(base.c))) return vc(c_pow_int(base.c, num));
        C out;
        if (den != 1 && c_pow_q(base.c, e, &out)) return vc(out);
        longjmp(need_series, 1);
    }
    if (den != 1) nm_fail("fractional powers of approximate numbers are not in this version");
    return vball(b_pow(as_ball(base), num, work_prec));
}

static Val arith(char op, Val a, Val b) {
    if (op == '^') return power_val(a, b);
    if (a.kind == V_Q && b.kind == V_Q) {
        switch (op) {
        case '+': return vq(q_add(a.q, b.q));
        case '-': return vq(q_sub(a.q, b.q));
        case '*': return vq(q_mul(a.q, b.q));
        case '/': return vq(q_div(a.q, b.q));
        }
    }
    if (is_exact(a) && is_exact(b)) {
        C x = as_c(a), y = as_c(b);
        switch (op) {
        case '+': return vc(c_add(x, y));
        case '-': return vc(c_sub(x, y));
        case '*': return vc(c_mul(x, y));
        case '/':
            if (!c_is_monomial(y) && c_has_plain(y)) longjmp(need_series, 1);
            return vc(c_div(x, y));                      /* by a single term, exactly: a a / x = a^2 x^-1 */
        }
    }
    Ball x = as_ball(a), y = as_ball(b);
    switch (op) {
    case '+': return vball(b_add(x, y, work_prec));
    case '-': return vball(b_sub(x, y, work_prec));
    case '*': return vball(b_mul(x, y, work_prec));
    case '/': return vball(b_div(x, y, work_prec));
    }
    nm_fail("internal: unknown operation");
    return a;
}

static Val name_val(const char *s, size_t len, int primes) {
    if (ovr.active && same(s, len, ovr.u, ovr.ul)) {
        if (primes > 7) nm_fail("derivatives above the seventh are not in this version");
        return vq(ovr.vals[primes]);
    }
    if (primes) nm_fail("%.*s' marks a derivative; it belongs inside an equation with a start", (int)len, s);
    Param *pp = param_find(s, len);
    if (pp) {
        if (pp->hasv) return pp->v;
        for (int i = 1; i < pp->d.a.n; i++)
            if (!c_is_zero(pp->d.a.c[i])) nm_fail("%.*s is a series here; a condition or an index needs a number", (int)len, s);
        if (pp->d.hasb) nm_fail("%.*s involves the unknown; a condition or an index needs a number", (int)len, s);
        return vc(pp->d.a.n ? pp->d.a.c[0] : c_zero());
    }
    Binding *b = lookup(s, len);
    if (b) {
        if (b->v.kind == V_FUNC) nm_fail("%s is a rule; write %s(...)", b->name, b->name);
        if (b->v.kind == V_SEQ) nm_fail("%s is a sequence; write %s[k]", b->name, b->name);
        if (b->v.kind == V_REC) longjmp(need_series, 1);
        if (b->v.kind == V_NUMREC) return eval(parse_text(b->v.src));   /* carried again to the places now asked */
        return b->v;
    }
    if (len == 1 && s[0] == 'i') return vc(c_imag_unit());   /* i: the square root of -1 */
    return vc(c_letter(letter_index(s, len)));       /* a letter */
}

static Val seq_term(Binding *b, int64_t k);
static int holds(Node *c);
static Val persist_val(Val v);

static Val call_rule(Binding *b, Node *n) {
    Def *d = b->v.def;
    if (n->nargs != d->np) nm_fail("%s takes %d argument%s", b->name, d->np, d->np == 1 ? "" : "s");
    Val *args = arena_alloc((size_t)d->np * sizeof(Val));
    for (int i = 0; i < d->np; i++) args[i] = eval(n->args[i]);
    if (++calldepth > 4000) nm_fail("%s calls itself more than 4000 deep; write it as a sequence, whose table is built forward", b->name);
    int saved_base = frame_base, saved_depth = pdepth;
    frame_base = pdepth;
    for (int i = 0; i < d->np; i++) { Param *p = param_push(d->pn[i], strlen(d->pn[i])); p->hasv = 1; p->v = args[i]; }
    Val v = eval(parse_text(d->src));
    pdepth = saved_depth; frame_base = saved_base; calldepth--;
    return v;
}

static int64_t whole(Val v, const char *what) {
    int64_t k;
    if (v.kind != V_Q || !q_is_int(v.q) || !z_fits_i64(v.q.num, &k)) nm_fail("%s must be a whole number", what);
    return k;
}

static Val eval(Node *n) {
    switch (n->kind) {
    case N_INDEX: {
        Binding *b = lookup(n->s, n->len);
        if (!b || b->v.kind != V_SEQ) nm_fail("%.*s is not a sequence (define it with let %.*s[0] = ..., %.*s[k] = ...)",
                                              (int)n->len, n->s, (int)n->len, n->s, (int)n->len, n->s);
        return seq_term(b, whole(eval(n->a), "an index"));
    }
    case N_SUM: {
        int64_t from = whole(eval(n->b), "the start of a sum"), to = whole(eval(n->c), "the end of a sum");
        if (to - from > 1000000) nm_fail("too many terms in a sum");
        Val acc = vq(qi(0));
        for (int64_t k = from; k <= to; k++) {
            Param *p = param_push(n->s, n->len);
            p->hasv = 1; p->v = vq(qi(k));
            Val t = eval(n->a);
            pdepth--;
            acc = arith('+', acc, t);
        }
        return acc;
    }
    case N_CASES:
        for (int i = 0; i < n->ncases; i++) if (holds(n->ccond[i])) return eval(n->cval[i]);
        return eval(n->c);
    case N_CMP: nm_fail("a comparison belongs in a definition by cases");
    case N_NUM: return vq(parse_number(n->s, n->len));
    case N_NAME: return name_val(n->s, n->len, n->primes);
    case N_NEG: return arith('-', vq(qi(0)), eval(n->a));
    case N_BIN: return arith(n->op, eval(n->a), eval(n->b));
    case N_SQRT: return power_val(eval(n->a), vq(q_make(z_from_i64(1), z_from_i64(2))));
    case N_APPLY: {
        Binding *b = lookup(n->s, n->len);
        if (b && b->v.kind == V_FUNC && !param_find(n->s, n->len)) return call_rule(b, n);
        if (n->nargs > 1) nm_fail("%.*s is not a rule of several arguments", (int)n->len, n->s);
        if (b && b->v.kind == V_REC) {
            Val arg = eval(n->a);
            if (arg.kind == V_POLY && c_has_plain(arg.c)) longjmp(need_series, 1);
            if (arg.kind == V_POLY && c_has_plain(arg.c)) longjmp(need_series, 1);
            return series_value(b, arg);
        }
        if (b && b->v.kind == V_POLY) {                     /* a polynomial in one letter, at a number */
            Val arg = eval(n->a);
            int li = -1, nl = 0;
            for (int l = 0; l < letter_count(); l++) if (c_uses(b->v.c, l)) { li = l; nl++; }
            if (arg.kind == V_Q && nl == 1) {
                Q acc = qi(0);
                for (int t = 0; t < b->v.c.nt; t++) {
                    CT term = b->v.c.t[t];
                    int64_t e;
                    if (!q_is_int(term.e[li]) || !z_fits_i64(term.e[li].num, &e)) nm_fail("a fractional power at a number is not in this version");
                    acc = q_add(acc, q_mul(term.k, q_pow(arg.q, e)));
                }
                return vq(acc);
            }
        }
        return arith('*', name_val(n->s, n->len, n->primes), eval(n->a));
    }
    case N_DERIV: case N_INTEG: {
        Val v = eval(n->a);
        if (!is_exact(v)) { if (n->kind == N_DERIV) return vq(qi(0)); longjmp(need_series, 1); }
        C p = as_c(v);
        int li = -1;
        if (n->len) li = letter_index(n->s, n->len);
        else {
            int nl = 0;
            for (int l = 0; l < letter_count(); l++) if (c_uses(p, l)) { li = l; nl++; }
            if (nl != 1) longjmp(need_series, 1);             /* the letter is decided by the series pass */
        }
        C r = c_zero();
        for (int t = 0; t < p.nt; t++) {
            CT term = p.t[t];
            Q e = term.e[li];
            if (n->kind == N_DERIV) {
                if (q_sign(e) == 0) continue;
                term.k = q_mul(term.k, e); term.e[li] = q_sub(e, qi(1));
            } else {                                       /* Newton's first rule: a x^m gives a x^(m+1)/(m+1) */
                if (q_cmp(e, qi(-1)) == 0) nm_fail("the fluent of 1/%s is a logarithm: write it as a root of an equation", letter_name(li));
                term.e[li] = q_add(e, qi(1)); term.k = q_div(term.k, term.e[li]);
            }
            C one; one.nt = 1; one.t = arena_alloc(sizeof(CT)); one.t[0] = term;
            r = c_add(r, one);
        }
        return vc(r);
    }
    case N_ROOT: {
        if (n->para) nm_fail("a root found with the parallelogram stands on its own line in this version");
        if (!n->c) longjmp(need_series, 1);
        Val l = eval(n->a), r = eval(n->b), g = eval(n->c);
        if (!is_exact(l) || !is_exact(r)) nm_fail("the equation must have exact coefficients");
        if (g.kind != V_Q) nm_fail("the starting value after 'near' must be an exact number");
        Poly p = c_to_poly(c_sub(as_c(l), as_c(r)), NULL);
        int64_t d0 = z_digits(g.q.den) + 2;
        Val v; memset(&v, 0, sizeof v);
        v.kind = V_ROOT;
        v.root = root_new(p, g.q, d0);
        return v;
    }
    }
    nm_fail("internal: unknown node");
    return vq(qi(0));
}

/* ---------------- the series pass ---------------- */

typedef struct {
    const char *var;            /* the letter of the series */
    int n;                      /* coefficients wanted */
    const char *unk; size_t unklen;
    Dual *yd; int r;            /* the unknown and its derivatives y, y', ..., y^(r) */
} SCtx;

static Dual dconst(Ser a) { Dual d; d.a = a; d.b = a; d.hasb = 0; return d; }

static Dual dadd(Dual x, Dual y, int sign) {
    Dual d;
    d.a = sign > 0 ? s_add(x.a, y.a) : s_sub(x.a, y.a);
    d.hasb = x.hasb || y.hasb;
    if (d.hasb) {
        Ser xb = x.hasb ? x.b : s_const(c_zero(), x.a.n), yb = y.hasb ? y.b : s_const(c_zero(), y.a.n);
        d.b = sign > 0 ? s_add(xb, yb) : s_sub(xb, yb);
    }
    return d;
}

static Dual dmul(Dual x, Dual y) {
    Dual d;
    d.a = s_mul(x.a, y.a);
    d.hasb = x.hasb || y.hasb;
    if (d.hasb) {
        Ser t = s_const(c_zero(), d.a.n);
        if (x.hasb) t = s_add(t, s_mul(x.b, y.a));
        if (y.hasb) t = s_add(t, s_mul(x.a, y.b));
        d.b = t;
    }
    return d;
}

static Dual ddiv(Dual x, Dual y) {
    Dual d;
    d.a = s_div(x.a, y.a);
    d.hasb = x.hasb || y.hasb;
    if (d.hasb) {                                           /* (b - (x/y) b') / y */
        Ser t = x.hasb ? x.b : s_const(c_zero(), x.a.n);
        if (y.hasb) t = s_sub(t, s_mul(d.a, y.b));
        d.b = s_div(t, y.a);
    }
    return d;
}

static Dual dpow(Dual x, Q alpha) {
    Dual d;
    d.a = s_pow_q(x.a, alpha);
    d.hasb = x.hasb;
    if (x.hasb) {
        Ser dfa;                                            /* x^(alpha - 1) */
        if (q_is_int(alpha) && q_sign(alpha) > 0) dfa = s_pow_q(x.a, q_sub(alpha, qi(1)));
        else dfa = s_div(d.a, x.a);
        d.b = s_mul(s_scale(dfa, c_const(alpha)), x.b);
    }
    return d;
}

static Dual sev(Node *n, SCtx *cx);
static Ser resolve(Node *root, const char *var, int n);

static C const_c_of(Dual d, const char *what) {
    for (int i = 1; i < d.a.n; i++)
        if (!c_is_zero(d.a.c[i])) nm_fail("%s must not depend on the series letter", what);
    if (d.hasb) nm_fail("%s must not involve the unknown", what);
    return d.a.n ? d.a.c[0] : c_zero();
}

static Q const_of(Dual d, const char *what) {
    Q k;
    if (!c_const_value(const_c_of(d, what), &k)) nm_fail("%s must be a number", what);
    return k;
}

/* a quantity in letters as a series in the letter var (whole non-negative powers) */
static Ser s_from_c(C c, const char *var, int n) {
    int vi = letter_index(var, strlen(var));
    Ser s = s_const(c_zero(), n);
    for (int t = 0; t < c.nt; t++) {
        CT term = c.t[t];
        int64_t e;
        if (!q_is_int(term.e[vi]) || !z_fits_i64(term.e[vi].num, &e) || e < 0)
            nm_fail("a negative or fractional power of %s: use the parallelogram (root of ... for y)", var);
        if (e >= n) continue;
        term.e[vi] = qi(0);
        C one; one.nt = 1; one.t = arena_alloc(sizeof(CT)); one.t[0] = term;
        s.c[e] = c_add(s.c[e], one);
    }
    return s;
}

static Ser recipe_series(Binding *b, const char *var, int n) {
    if (strcmp(b->v.var, var)) nm_fail("%s is a series in %s; write %s(%s) to substitute", b->name, b->v.var, b->name, var);
    SCtx cx; memset(&cx, 0, sizeof cx);
    cx.var = b->v.var; cx.n = n;
    return sev(parse_text(b->v.src), &cx).a;
}

static Dual sname(const char *s, size_t len, int primes, SCtx *cx) {
    if (cx->unk && same(s, len, cx->unk, cx->unklen)) {
        if (primes > cx->r) nm_fail("internal: derivative order");
        return cx->yd[primes];
    }
    if (primes) nm_fail("%.*s' marks a derivative of an unknown with no start", (int)len, s);
    Param *pp = param_find(s, len);
    if (pp) {
        if (pp->hasd) return pp->d;
        if (pp->v.kind == V_Q) return dconst(s_const(c_const(pp->v.q), cx->n));
        if (pp->v.kind == V_POLY) return dconst(s_from_c(pp->v.c, cx->var, cx->n));
        nm_fail("%.*s is an approximate number; it cannot enter a series", (int)len, s);
    }
    if (same(s, len, cx->var, strlen(cx->var))) return dconst(s_var(cx->n));
    Binding *b = lookup(s, len);
    if (!b && len == 1 && s[0] == 'i') return dconst(s_const(c_imag_unit(), cx->n));
    if (!b) return dconst(s_const(c_letter(letter_index(s, len)), cx->n));    /* a given letter: a, b, ... */
    switch (b->v.kind) {
    case V_Q: return dconst(s_const(c_const(b->v.q), cx->n));
    case V_POLY: return dconst(s_from_c(b->v.c, cx->var, cx->n));
    case V_REC: return dconst(recipe_series(b, cx->var, cx->n));
    default: nm_fail("%s is an irrational number; irrational coefficients come later", b->name);
    }
    return dconst(s_const(c_zero(), cx->n));
}

static Dual sev(Node *n, SCtx *cx) {
    switch (n->kind) {
    case N_INDEX: {
        Val v = eval(n);
        if (v.kind == V_Q) return dconst(s_const(c_const(v.q), cx->n));
        return dconst(s_from_c(v.c, cx->var, cx->n));
    }
    case N_SUM: {
        int64_t from = whole(eval(n->b), "the start of a sum"), to = whole(eval(n->c), "the end of a sum");
        if (to - from > 1000000) nm_fail("too many terms in a sum");
        Dual acc = dconst(s_const(c_zero(), cx->n));
        for (int64_t k = from; k <= to; k++) {
            Param *p = param_push(n->s, n->len);
            p->hasv = 1; p->v = vq(qi(k));
            Dual t = sev(n->a, cx);
            pdepth--;
            acc = dadd(acc, t, 1);
        }
        return acc;
    }
    case N_CASES:
        for (int i = 0; i < n->ncases; i++) if (holds(n->ccond[i])) return sev(n->cval[i], cx);
        return sev(n->c, cx);
    case N_CMP: nm_fail("a comparison belongs in a definition by cases");
    case N_NUM: return dconst(s_const(c_const(parse_number(n->s, n->len)), cx->n));
    case N_NAME: return sname(n->s, n->len, n->primes, cx);
    case N_NEG: return dadd(dconst(s_const(c_zero(), cx->n)), sev(n->a, cx), -1);
    case N_BIN: {
        Dual x = sev(n->a, cx);
        if (n->op == '^') return dpow(x, const_of(sev(n->b, cx), "an exponent"));
        Dual y = sev(n->b, cx);
        switch (n->op) {
        case '+': return dadd(x, y, 1);
        case '-': return dadd(x, y, -1);
        case '*': return dmul(x, y);
        case '/': return ddiv(x, y);
        }
        break;
    }
    case N_SQRT: return dpow(sev(n->a, cx), q_make(z_from_i64(1), z_from_i64(2)));
    case N_APPLY: {
        Binding *b = lookup(n->s, n->len);
        if (b && b->v.kind == V_FUNC && !param_find(n->s, n->len)) {
            Def *d = b->v.def;
            if (n->nargs != d->np) nm_fail("%s takes %d argument%s", b->name, d->np, d->np == 1 ? "" : "s");
            Dual *args = arena_alloc((size_t)d->np * sizeof(Dual));
            for (int i = 0; i < d->np; i++) args[i] = sev(n->args[i], cx);
            if (++calldepth > 4000) nm_fail("%s calls itself more than 4000 deep; write it as a sequence, whose table is built forward", b->name);
            int saved_base = frame_base, saved_depth = pdepth;
            frame_base = pdepth;
            for (int i = 0; i < d->np; i++) { Param *p = param_push(d->pn[i], strlen(d->pn[i])); p->hasd = 1; p->d = args[i]; }
            Dual r = sev(parse_text(d->src), cx);
            pdepth = saved_depth; frame_base = saved_base; calldepth--;
            return r;
        }
        if (n->nargs > 1) nm_fail("%.*s is not a rule of several arguments", (int)n->len, n->s);
        Dual g = sev(n->a, cx);
        if (!b || b->v.kind != V_REC) return dmul(sname(n->s, n->len, n->primes, cx), g);
        Ser f = recipe_series(b, b->v.var, g.a.n);                /* f(g), and its moment f'(g) g' */
        Dual d;
        d.a = s_compose(f, g.a);
        d.hasb = g.hasb;
        if (g.hasb) d.b = s_mul(s_compose(s_deriv(f), s_trunc(g.a, f.n - 1)), g.b);
        return d;
    }
    case N_DERIV: case N_INTEG: {
        if (n->len && !same(n->s, n->len, cx->var, strlen(cx->var))) {
            if (n->kind == N_DERIV) return dconst(s_const(c_zero(), cx->n));
            nm_fail("integral in %.*s inside a series in %s", (int)n->len, n->s, cx->var);
        }
        Dual x = sev(n->a, cx), d;
        d.hasb = x.hasb;
        if (n->kind == N_DERIV) { d.a = s_deriv(x.a); d.b = x.hasb ? s_deriv(x.b) : d.a; }
        else {
            d.a = s_trunc(s_integ(x.a), cx->n);
            d.b = x.hasb ? s_trunc(s_integ(x.b), cx->n) : d.a;
        }
        return d;
    }
    case N_ROOT:
        if (n->para) nm_fail("a root found with the parallelogram stands on its own line in this version");
        if (n->c) nm_fail("a number found with 'near' cannot be a series coefficient unless it is rational");
        return dconst(resolve(n, cx->var, cx->n));
    }
    nm_fail("internal: unknown node");
    return dconst(s_const(c_zero(), 0));
}

static int max_primes(Node *n, const char *u, size_t ul) {
    if (!n) return 0;
    int m = 0;
    if ((n->kind == N_NAME || n->kind == N_APPLY) && same(n->s, n->len, u, ul)) m = n->primes;
    int a = max_primes(n->a, u, ul), b = max_primes(n->b, u, ul), c = max_primes(n->c, u, ul);
    if (a > m) m = a;
    if (b > m) m = b;
    if (c > m) m = c;
    return m;
}

/* Resolution of an equation for a series y, term by term from its start (Methodus, Problem 2; De analysi).
 * Put y = (terms found) + o x^d; the lowest term of the equation that o reaches fixes the new coefficient.
 * At the end the whole series is substituted back, and the equation must vanish to the order claimed. */
static Ser resolve(Node *root, const char *var, int n) {
    Cond *cs = root->conds;
    const char *u = cs[0].s; size_t ul = cs[0].len;
    for (int i = 1; i < root->nconds; i++)
        if (!same(cs[i].s, cs[i].len, u, ul)) nm_fail("all starts must be for the same unknown");
    int r = max_primes(root->a, u, ul);
    int r2 = max_primes(root->b, u, ul);
    if (r2 > r) r = r2;
    if (n > 2000) nm_fail("order too large");
    /* the start: y(0), y'(0), ..., y^(r-1)(0); for an equation without derivatives, y(0) chooses the branch */
    int need = r ? r : 1;
    C *c = arena_alloc((size_t)(n + r + 1) * sizeof(C));
    for (int i = 0; i < n + r + 1; i++) c[i] = c_zero();
    int *given = arena_alloc((size_t)need * sizeof(int));
    memset(given, 0, (size_t)need * sizeof(int));
    SCtx k0; memset(&k0, 0, sizeof k0); k0.var = var; k0.n = 1;
    for (int i = 0; i < root->nconds; i++) {
        if (q_sign(const_of(sev(cs[i].at, &k0), "the point of a start")) != 0)
            nm_fail("starts must be given at 0 in this version");
        int j = cs[i].primes;
        if (j >= need) nm_fail("the equation has order %d; %.*s with %d marks is not a start for it", r, (int)ul, u, j);
        C v = const_c_of(sev(cs[i].val, &k0), "a starting value");
        Q fact = qi(1);
        for (int f = 2; f <= j; f++) fact = q_mul(fact, qi(f));
        c[j] = c_scale(v, q_div(qi(1), fact));
        given[j] = 1;
    }
    for (int j = 0; j < need; j++)
        if (!given[j]) nm_fail("missing start: %.*s%.*s(0)", (int)ul, u, j, "''''''''''''''''");
    Node *eqn = bin('-', root->a, root->b);
    SCtx cx; memset(&cx, 0, sizeof cx);
    cx.var = var; cx.unk = u; cx.unklen = ul; cx.r = r;
    cx.yd = arena_alloc((size_t)(r + 1) * sizeof(Dual));
    for (int d = need; d < n; d++) {
        int k = d - r;                                  /* the equation's term that decides y's term of degree d */
        int m = k + 1 + r;
        Ser Y; Y.n = m; Y.c = arena_alloc((size_t)m * sizeof(C));
        for (int i = 0; i < m; i++) Y.c[i] = i < d ? c[i] : c_zero();
        Ser O = s_const(c_zero(), m);
        O.c[d] = c_const(qi(1));
        for (int j = 0; j <= r; j++) {
            cx.yd[j].a = s_trunc(Y, k + 1); cx.yd[j].b = s_trunc(O, k + 1); cx.yd[j].hasb = 1;
            Y = s_deriv(Y); O = s_deriv(O);
        }
        cx.n = k + 1;
        Dual e = sev(eqn, &cx);
        C F = e.a.n > k ? e.a.c[k] : c_zero();
        C G = (e.hasb && e.b.n > k) ? e.b.c[k] : c_zero();
        if (c_is_zero(G))
            nm_fail("the start is a multiple root: the next term is not fixed by the lowest terms (use the parallelogram: root of ... for %.*s)", (int)ul, u);
        c[d] = c_neg(c_div(F, G));
    }
    /* substitute back */
    Ser Y; Y.n = n; Y.c = c;
    int check = n - r;
    for (int j = 0; j <= r; j++) { cx.yd[j] = dconst(s_trunc(Y, check)); Y = s_deriv(Y); }
    cx.n = check;
    Dual e = sev(eqn, &cx);
    for (int i = 0; i < e.a.n && i < check; i++)
        if (!c_is_zero(e.a.c[i])) {
            if (i == 0 && r == 0) nm_fail("%.*s(0) = %s does not satisfy the equation at %s = 0", (int)ul, u, c_to_str(c[0]), var);
            nm_fail("substituting the series back leaves a remainder at %s^%d: the equation cannot be resolved term by term here", var, i);
        }
    Ser out; out.n = n; out.c = c;
    return out;
}

/* the free letters of a statement (not the unknown of an equation); the series letter is one of them */
typedef struct { const char *s[NM_MAXL + 4]; size_t len[NM_MAXL + 4]; int n; } Letters;

static void add_letter(Letters *L, const char *s, size_t len) {
    for (int i = 0; i < L->n; i++) if (same(L->s[i], L->len[i], s, len)) return;
    if (L->n < NM_MAXL + 4) { L->s[L->n] = s; L->len[L->n] = len; L->n++; }
}

static const char *sum_skip[64]; static size_t sum_skiplen[64]; static int nsum_skip;

static void letters(Node *n, Letters *L, const char *skip, size_t skiplen) {
    if (!n) return;
    if (n->kind == N_SUM) {
        letters(n->b, L, skip, skiplen); letters(n->c, L, skip, skiplen);
        if (nsum_skip < 64) { sum_skip[nsum_skip] = n->s; sum_skiplen[nsum_skip] = n->len; nsum_skip++; }
        letters(n->a, L, skip, skiplen);
        nsum_skip--;
        return;
    }
    if (n->kind == N_CASES) {
        for (int i = 0; i < n->ncases; i++) letters(n->cval[i], L, skip, skiplen);
        letters(n->c, L, skip, skiplen);
        return;
    }
    if (n->kind == N_NAME)
        for (int i = 0; i < nsum_skip; i++) if (same(n->s, n->len, sum_skip[i], sum_skiplen[i])) return;
    if (n->kind == N_INDEX) { letters(n->a, L, skip, skiplen); return; }
    if (n->kind == N_APPLY && n->nargs > 1) {
        for (int i = 0; i < n->nargs; i++) letters(n->args[i], L, skip, skiplen);
        return;
    }
    if (n->kind == N_NAME || n->kind == N_APPLY) {
        Binding *b = lookup(n->s, n->len);
        if (b && (b->v.kind == V_FUNC || b->v.kind == V_SEQ)) { letters(n->a, L, skip, skiplen); return; }
        if (b && b->v.kind == V_REC && n->kind == N_NAME) add_letter(L, b->v.var, strlen(b->v.var));
        else if (b && b->v.kind == V_POLY) {
            for (int l = 0; l < letter_count(); l++) if (c_uses(b->v.c, l)) add_letter(L, letter_name(l), strlen(letter_name(l)));
        } else if (!b && !(skip && same(n->s, n->len, skip, skiplen)) && !(n->len == 1 && n->s[0] == 'i'))
            add_letter(L, n->s, n->len);
    }
    if ((n->kind == N_DERIV || n->kind == N_INTEG) && n->len) add_letter(L, n->s, n->len);
    if (n->kind == N_ROOT && n->nconds) { skip = n->conds[0].s; skiplen = n->conds[0].len; }
    if (n->kind == N_ROOT && n->para) { skip = n->us; skiplen = n->ul; }
    letters(n->a, L, skip, skiplen);
    letters(n->b, L, skip, skiplen);
    letters(n->c, L, skip, skiplen);
    letters(n->start, L, skip, skiplen);
}

/* ---------------- the value of a series at a number ---------------- */

/* The term rule. If the equation is  sum_j p_j(x) y^(j) + h(x) = 0  with polynomials p_j and h, then putting
 * y = sum c_k x^k gives, for n large enough,
 *     D(n) c_n = sum_{t>=1} N_t(n) c_{n-t} - h_{n-s},
 * each term from the ones before, as Newton's A, B, C, D. */
typedef struct { int s, T; Poly D; Poly *N; Poly h; } Rule;

static Poly poly_of(Val v, const char *var) {
    if (v.kind == V_Q) return p_const(v.q);
    if (v.kind == V_POLY) { Poly p = c_to_poly(v.c, var); p.var = var; return p; }
    nm_fail("the equation's coefficients must be polynomials in %s", var);
    return p_const(qi(0));
}

static Poly falling(Poly n, int64_t shift, int j) {   /* (n - shift)(n - shift - 1) ... j factors */
    Poly r = p_const(qi(1));
    for (int u = 0; u < j; u++) r = p_mul(r, p_sub(n, p_const(qi(shift + u))));
    return r;
}

static Q p_at(Poly p, Q x) {
    Q acc = qi(0);
    for (int i = p.deg; i >= 0; i--) acc = q_add(q_mul(acc, x), p.c[i]);
    return acc;
}

static int p_equal(Poly a, Poly b) { return p_sub(a, b).deg < 0; }

static Rule term_rule(Node *root, const char *var, int r) {
    Node *E = bin('-', root->a, root->b);
    jmp_buf saved;
    memcpy(saved, need_series, sizeof saved);
    if (setjmp(need_series)) {
        memcpy(need_series, saved, sizeof saved);
        ovr.active = 0;
        nm_fail("a value needs the equation with polynomial coefficients (multiply out the denominators)");
    }
    ovr.active = 1; ovr.u = root->conds[0].s; ovr.ul = root->conds[0].len;
    for (int j = 0; j < 8; j++) ovr.vals[j] = qi(0);
    Poly h = poly_of(eval(E), var);
    Poly *p = arena_alloc((size_t)(r + 1) * sizeof(Poly));
    for (int j = 0; j <= r; j++) {                           /* read off each coefficient, and check linearity */
        ovr.vals[j] = qi(1);
        p[j] = p_sub(poly_of(eval(E), var), h);
        ovr.vals[j] = qi(2);
        Poly twice = p_sub(poly_of(eval(E), var), h);
        ovr.vals[j] = qi(0);
        if (!p_equal(twice, p_scale(p[j], qi(2)))) goto nonlinear;
        for (int i = 0; i < j; i++) {
            ovr.vals[i] = qi(1); ovr.vals[j] = qi(1);
            Poly both = p_sub(poly_of(eval(E), var), h);
            ovr.vals[i] = qi(0); ovr.vals[j] = qi(0);
            if (!p_equal(both, p_add(p[i], p[j]))) goto nonlinear;
        }
    }
    ovr.active = 0;
    memcpy(need_series, saved, sizeof saved);
    Rule R; memset(&R, 0, sizeof R);
    R.s = -1000000;
    int lowshift = 1000000;
    for (int j = 0; j <= r; j++)
        for (int i = 0; i <= p[j].deg; i++)
            if (q_sign(p[j].c[i])) { if (j - i > R.s) R.s = j - i; if (j - i < lowshift) lowshift = j - i; }
    if (R.s == -1000000) nm_fail("the equation does not involve the unknown");
    R.T = R.s - lowshift;
    Poly nvar = p_var("n");
    R.D = p_const(qi(0));
    R.N = arena_alloc((size_t)(R.T + 1) * sizeof(Poly));
    for (int t = 0; t <= R.T; t++) R.N[t] = p_const(qi(0));
    for (int j = 0; j <= r; j++)
        for (int i = 0; i <= p[j].deg; i++) {
            if (!q_sign(p[j].c[i])) continue;
            int t = R.s - (j - i);
            Poly term = p_scale(falling(nvar, t, j), p[j].c[i]);
            if (t == 0) R.D = p_add(R.D, term);
            else R.N[t] = p_sub(R.N[t], term);
        }
    R.h = h;
    return R;
nonlinear:
    ovr.active = 0;
    memcpy(need_series, saved, sizeof saved);
    nm_fail("a value needs an equation that is linear in the unknown and its derivatives");
    return (Rule){0};
}

static Q qabs(Q a) { return q_sign(a) < 0 ? q_neg(a) : a; }

/* an upper bound for |P(n)/D(n)| valid for every n >= N, when deg P <= deg D; negative if none yet */
static Q ratio_bound(Poly P, Poly D, int64_t N) {
    int d = D.deg;
    Q nn = qi(N), num = qi(0), den = qabs(D.c[d]);
    for (int i = 0; i <= P.deg; i++) num = q_add(num, q_div(qabs(P.c[i]), q_pow(nn, d - i)));
    for (int i = 0; i < d; i++) den = q_sub(den, q_div(qabs(D.c[i]), q_pow(nn, d - i)));
    if (q_sign(den) <= 0) return qi(-1);
    return q_div(num, den);
}

static Ball ball_widen(Ball b, Q err) {               /* add |err| to the radius */
    Z num = z_abs(err.num), den = err.den, q, rem;
    if (b.e < 0) num = z_mul_pow10(num, -b.e); else den = z_mul(den, z_pow10(b.e));
    z_divmod(num, den, &q, &rem);
    b.r = z_add(z_add(b.r, q), z_from_i64(rem.s ? 1 : 0));
    return b;
}

static Val series_value(Binding *b, Val arg) {
    Q a = arg.kind == V_Q ? arg.q : q_from_z(z_zero());
    Node *root = parse_text(b->v.src);
    if (root->kind != N_ROOT || !root->nconds)
        nm_fail("the value of %s needs its definition as an equation with a start (root of ..., y(0) = ...)", b->name);
    const char *u = root->conds[0].s; size_t ul = root->conds[0].len;
    int r = max_primes(root->a, u, ul), r2 = max_primes(root->b, u, ul);
    if (r2 > r) r = r2;
    Rule R = term_rule(root, b->v.var, r);
    if (R.D.deg < 0) nm_fail("the equation fixes no term of %s", b->name);
    Q rr = qabs(a);
    Ball Aball;
    if (arg.kind != V_Q) {                               /* an approximate argument: the tail with an upper bound of |a| */
        Aball = as_ball(arg);
        Z top = z_add(z_abs(Aball.m), Aball.r);
        rr = Aball.e >= 0 ? q_from_z(z_mul_pow10(top, Aball.e)) : q_make(top, z_pow10(-Aball.e));
    }
    /* how fast the terms shrink in the end: sum over t of lim |N_t / D| r^t must be below 1 */
    Q ginf = qi(0);
    for (int t = 1; t <= R.T; t++) {
        if (R.N[t].deg > R.D.deg) nm_fail("the terms of %s grow too fast: the series has no value away from 0", b->name);
        if (R.N[t].deg == R.D.deg) ginf = q_add(ginf, q_mul(qabs(q_div(R.N[t].c[R.D.deg], R.D.c[R.D.deg])), q_pow(rr, t)));
    }
    if (q_cmp_one(ginf) >= 0)
        nm_fail("at %s the terms of %s shrink too slowly to bound the rest; use smaller arguments and exact relations, as Newton did with 0.1 and 0.2",
                arg.kind == V_Q ? q_to_str(a) : "this argument", b->name);
    /* from where on the rule holds: beyond the integer roots of D and the inhomogeneous terms */
    Q cb = qi(0);
    for (int i = 0; i < R.D.deg; i++) {
        Q v = qabs(q_div(R.D.c[i], R.D.c[R.D.deg]));
        if (q_cmp(v, cb) > 0) cb = v;
    }
    Z cbq, cbr;
    z_divmod(cb.num, cb.den, &cbq, &cbr);
    int64_t start;
    if (!z_fits_i64(cbq, &start) || start > 100000) nm_fail("the rule for %s starts too late", b->name);
    start += 2;
    if (start < R.s + R.h.deg + 1) start = R.s + R.h.deg + 1;
    if (start < R.T) start = R.T;
    /* seeds from the resolution; the rule must reproduce the next few, as a check */
    int nseed = (int)start + 4;
    Ser seed = resolve(root, b->v.var, nseed);
    int cap = nseed + 64;
    Q *c = arena_alloc((size_t)cap * sizeof(Q));
    for (int i = 0; i < nseed; i++)
        if (!c_const_value(seed.c[i], &c[i])) nm_fail("a value needs numbers; %s has letters in its terms", b->name);
    int ncoef = (int)start;
    for (int n = (int)start; n < nseed; n++) {
        Q acc = qi(0);
        for (int t = 1; t <= R.T && t <= n; t++) acc = q_add(acc, q_mul(p_at(R.N[t], qi(n)), c[n - t]));
        int m = n - R.s;
        if (m >= 0 && m <= R.h.deg) acc = q_sub(acc, R.h.c[m]);
        Q cn = q_div(acc, p_at(R.D, qi(n)));
        if (q_cmp(cn, c[n]) != 0) nm_fail("internal check failed: the term rule of %s disagrees with its resolution (vitiose)", b->name);
    }
    ncoef = nseed;
    if (arg.kind == V_Q && q_sign(a) == 0) return vq(c[0]);
    /* sum term by term in places (the coefficients exact, the powers of a carried as balls); stop when the rest
     * is certainly below the places asked */
    int64_t want = work_prec - GUARD_DIGITS + 10;
    Q eps = q_make(z_from_i64(1), z_pow10(want));
    Ball S = b_from_q(qi(0), work_prec), A = arg.kind == V_Q ? b_from_q(a, work_prec) : Aball, P = b_from_q(qi(1), work_prec);
    Ball *term = arena_alloc((size_t)cap * sizeof(Ball));
    for (int k = 0;; k++) {
        if (k >= ncoef) {
            if (ncoef == cap) {
                int ncap = cap * 2;
                Q *nc = arena_alloc((size_t)ncap * sizeof(Q));
                Ball *nt = arena_alloc((size_t)ncap * sizeof(Ball));
                memcpy(nc, c, (size_t)cap * sizeof(Q)); memcpy(nt, term, (size_t)cap * sizeof(Ball));
                c = nc; term = nt; cap = ncap;
            }
            int n = ncoef;
            Q acc = qi(0);
            for (int t = 1; t <= R.T; t++) acc = q_add(acc, q_mul(p_at(R.N[t], qi(n)), c[n - t]));
            int m = n - R.s;
            if (m >= 0 && m <= R.h.deg) acc = q_sub(acc, R.h.c[m]);
            c[n] = q_div(acc, p_at(R.D, qi(n)));
            ncoef++;
        }
        term[k] = q_sign(c[k]) ? b_mul(b_from_q(c[k], work_prec), P, work_prec) : b_from_q(qi(0), work_prec);
        S = b_add(S, term[k], work_prec);
        P = b_mul(P, A, work_prec);
        int64_t N = k + 1;
        if (N < start || N < R.T) continue;
        /* the careful bound only once the last terms are small */
        int small = 1;
        for (int t = 1; t <= R.T; t++) {
            Ball b = term[N - t];
            Z top = z_add(z_abs(b.m), b.r);
            if (top.s && b.e + z_digits(top) > -want + 2) { small = 0; break; }
        }
        if (!small) { if (k > 2000000) nm_fail("too many terms"); continue; }
        Q gamma = qi(0);
        int ok = 1;
        for (int t = 1; t <= R.T; t++) {
            Q u = ratio_bound(R.N[t], R.D, N);
            if (q_sign(u) < 0) { ok = 0; break; }
            gamma = q_add(gamma, q_mul(u, q_pow(rr, t)));
        }
        if (!ok || q_cmp_one(gamma) >= 0) continue;
        Q W = qi(0);
        for (int t = 1; t <= R.T; t++) {                 /* an upper bound of |term|: (|m| + r) 10^e */
            Ball b = term[N - t];
            Z top = z_add(z_abs(b.m), b.r);
            Q up = b.e >= 0 ? q_from_z(z_mul_pow10(top, b.e)) : q_make(top, z_pow10(-b.e));
            if (q_cmp(up, W) > 0) W = up;
        }
        Q tail = q_div(q_mul(q_mul(qi(R.T), gamma), W), q_sub(qi(1), gamma));
        if (q_cmp(tail, eps) <= 0) return vball(ball_widen(S, tail));
    }
}

/* ---------------- conditions and sequences ---------------- */

static int holds(Node *c) {
    Val a = eval(c->a), b = eval(c->b);
    if (a.kind == V_Q && b.kind == V_Q) {
        int k = q_cmp(a.q, b.q);
        switch (c->op) {
        case '=': return k == 0;
        case '<': return k < 0;
        case '>': return k > 0;
        case 'l': return k <= 0;
        case 'g': return k >= 0;
        }
    }
    if (c->op == '=' && is_exact(a) && is_exact(b)) return c_equal(as_c(a), as_c(b));
    nm_fail("a condition compares exact numbers (a letter or an approximate number cannot be ordered here)");
    return 0;
}

/* the k-th term of a sequence, from a table built forward, as Newton built his tables */
static Val seq_term(Binding *b, int64_t k) {
    Def *d = b->v.def;
    if (k < 0) nm_fail("%s[%lld]: indices start at 0", b->name, (long long)k);
    if (d->gen != let_gen) { d->gen = let_gen; d->known = 0; d->building = 0; }
    if (k < d->known) return d->table[k];
    if (d->building == stmt_id && k >= d->progress)
        nm_fail("%s[%lld] is needed to find itself: a term must come from earlier terms", b->name, (long long)k);
    if (k > 10000000) nm_fail("index too large");
    long saved_building = d->building; int64_t saved_progress = d->progress;
    for (int64_t i = d->known; i <= k; i++) {
        d->building = stmt_id; d->progress = i;
        Val v;
        int init = -1;
        for (int j = 0; j < d->ninit; j++) if (d->iidx[j] == i) init = j;
        int saved_base = frame_base, saved_depth = pdepth;
        frame_base = pdepth;
        if (init >= 0) v = eval(parse_text(d->isrc[init]));
        else {
            int64_t first = 0;
            for (int j = 0; j < d->ninit; j++) if (d->iidx[j] + 1 > first) first = d->iidx[j] + 1;
            if (i < first) nm_fail("%s[%lld] has no starting value", b->name, (long long)i);
            Param *p = param_push(d->pn[0], strlen(d->pn[0]));
            p->hasv = 1; p->v = vq(qi(i));
            v = eval(parse_text(d->src));
        }
        pdepth = saved_depth; frame_base = saved_base;
        if (!is_exact(v)) nm_fail("the terms of %s must be exact (rational numbers or letters)", b->name);
        if (d->known == d->cap) {
            int ncap = d->cap ? d->cap * 2 : 64;
            Val *nt = perm_alloc((size_t)ncap * sizeof(Val));
            if (d->known) memcpy(nt, d->table, (size_t)d->known * sizeof(Val));
            d->table = nt; d->cap = ncap;
        }
        d->table[d->known++] = persist_val(v);
    }
    d->building = saved_building; d->progress = saved_progress;
    return d->table[k];
}

/* ---------------- writing results ---------------- */

static char *show(Val v, int64_t places, int asked);
static char *show(Val v, int64_t places, int asked) {
    char *out = arena_alloc(4096), *body;
    switch (v.kind) {
    case V_Q:
        if (!asked) { body = q_to_str(v.q); snprintf(out, 4096, "%s", "[exact]"); break; }
        {
            Z s = z_div_round(z_mul_pow10(v.q.num, places), v.q.den);
            body = fixed_str(s, places);
            Z q, rem;
            z_divmod(z_mul_pow10(v.q.num, places), v.q.den, &q, &rem);
            if (rem.s == 0) snprintf(out, 4096, "[exact]");
            else if (strlen(q_to_str(v.q)) > 60) snprintf(out, 4096, "[exact fraction, rounded to %lld places]", (long long)places);
            else snprintf(out, 4096, "[exact value %s, rounded to %lld places]", q_to_str(v.q), (long long)places);
        }
        break;
    case V_POLY: {
        C re, im;
        int lone = v.c.nt == 1 && q_is_int(v.c.t[0].k) && z_is_one(v.c.t[0].k.num);
        int li = -1, nl = 0;
        for (int l = 0; l < letter_count(); l++) if (c_uses(v.c, l)) { li = l; nl++; }
        if (asked && lone && nl == 1 && letter_is_surd(li) && surd_root(li) && q_cmp_one(v.c.t[0].e[li]) == 0) {
            Val rv; memset(&rv, 0, sizeof rv); rv.kind = V_ROOT; rv.root = surd_root(li);
            return show(rv, places, asked);                /* a named root: its certified places */
        }
        if (!asked || c_has_plain(v.c)) { body = c_to_str(v.c); snprintf(out, 4096, "[exact]"); break; }
        if (split_imag(v.c, &re, &im)) {                   /* a + b i */
            Ball br = ball_of_c(re, work_prec), bi = ball_of_c(im, work_prec);
            int64_t kr = b_guaranteed_places(br, places), ki = b_guaranteed_places(bi, places);
            int64_t k = kr < ki ? kr : ki;
            if (k < 0) nm_fail("not even the units place is guaranteed");
            char *sr = fixed_str(b_round_to_places(br, k), k), *si = fixed_str(b_round_to_places(bi, k), k);
            body = arena_alloc(strlen(sr) + strlen(si) + 8);
            if (si[0] == '-') sprintf(body, "%s - %si", sr, si + 1); else sprintf(body, "%s + %si", sr, si);
            if (k == places) snprintf(out, 4096, "[bounded: %lld places guaranteed]", (long long)k);
            else snprintf(out, 4096, "[bounded: only %lld of %lld places guaranteed]", (long long)k, (long long)places);
            break;
        }
        return show(vball(ball_of_c(v.c, work_prec)), places, asked);
    }
    case V_ROOT: {
        Root *r = v.root;
        root_refine(r, places);
        body = fixed_str(z_div_round(r->X, z_pow10(r->D - places)), places);
        if (!r->certified)
            snprintf(out, 4096, "[not certified: %s shows no sign change here; a double root?]", root_equation_str(r));
        else
            snprintf(out, 4096, "[certified: sign change of %s, %lld places]", root_equation_str(r), (long long)places);
        break;
    }
    case V_BALL: {
        int64_t k = b_guaranteed_places(v.ball, places);
        if (k < 0) {
            body = fixed_str(b_round_to_places(v.ball, 0), 0);
            snprintf(out, 4096, "[bounded: not even the units place is guaranteed]");
        } else {
            body = fixed_str(b_round_to_places(v.ball, k), k);
            if (k == places) snprintf(out, 4096, "[bounded: %lld places guaranteed]", (long long)k);
            else snprintf(out, 4096, "[bounded: only %lld of %lld places guaranteed]", (long long)k, (long long)places);
        }
        break;
    }
    default: body = "?";
    }
    char *line = arena_alloc(strlen(body) + strlen(out) + 4);
    sprintf(line, "%s  %s", body, out);
    return line;
}

static Val persist_val(Val v) {
    switch (v.kind) {
    case V_Q: v.q = q_persist(v.q); break;
    case V_POLY: v.c = c_persist(v.c); break;
    case V_BALL: v.ball.m = z_persist(v.ball.m); v.ball.r = z_persist(v.ball.r); break;
    default: break;                                     /* roots and recipes already live in permanent memory */
    }
    return v;
}

static void bind(const char *name, size_t len, Val v) {
    let_gen++;
    Binding *b = perm_alloc(sizeof *b);
    b->name = dup_perm(name, len);
    b->v = persist_val(v);
    b->next = names; names = b;
}

/* root of F = 0 for y [starting y = S] to x^N */
static char *para_root(Node *e, const char *to_var, size_t to_len, int64_t order) {
    if (!e->us) nm_fail("name the unknown: root of ... for y");
    if (max_primes(e->a, e->us, e->ul) || max_primes(e->b, e->us, e->ul))
        nm_fail("the parallelogram is for equations without derivatives; give a start y(0) = ... for a fluxional one");
    Val l, r;
    jmp_buf saved;
    memcpy(saved, need_series, sizeof saved);
    if (setjmp(need_series)) {
        memcpy(need_series, saved, sizeof saved);
        nm_fail("the parallelogram needs a polynomial equation (multiply out the denominators)");
    }
    l = eval(e->a); r = eval(e->b);
    C F = c_sub(as_c(l), as_c(r));
    C S = c_zero();
    int have = 0;
    if (e->start) { S = as_c(eval(e->start)); have = 1; }
    memcpy(need_series, saved, sizeof saved);
    int yi = letter_index(e->us, e->ul);
    int xi = -1;
    if (to_var) xi = letter_index(to_var, to_len);
    else {
        int n = 0, nf = 0, xf = -1;
        for (int li = 0; li < letter_count(); li++)
            if (li != yi && letter_name(li)[0] != '#' && c_uses(F, li)) {
                xi = li; n++;
                if (flowing_letter(letter_name(li), strlen(letter_name(li)))) { xf = li; nf++; }
            }
        if (n > 1 && nf == 1) xi = xf;                 /* Newton: a, b, c are given; x, z flow */
        else if (n != 1) nm_fail("name the letter of the series with 'to x^8'");
    }
    if (c_uses(S, yi)) nm_fail("the start must not contain %.*s", (int)e->ul, e->us);
    return parallelogram(F, xi, yi, &S, have, order);
}

#ifndef NM_LIBDIR
#define NM_LIBDIR "lib"
#endif

char *nm_run(const char *line, int *failed);

static char *use_library(const char *name, size_t len) {
    const char *dir = getenv("NEWTONMATH_LIB");
    char path[1024];
    snprintf(path, sizeof path, "%s/%.*s.nm", dir ? dir : NM_LIBDIR, (int)len, name);
    FILE *f = fopen(path, "r");
    if (!f) nm_fail("cannot open %s", path);
    jmp_buf saved;
    memcpy(saved, nm_on_error, sizeof saved);
    char *buf = NULL; size_t cap = 0;
    int defs = 0, lineno = 0;
    char *err = NULL;
    while (getline(&buf, &cap, f) >= 0) {
        lineno++;
        int bad;
        char *out = nm_run(buf, &bad);
        if (bad) {
            err = arena_alloc(strlen(out) + strlen(path) + 32);
            sprintf(err, "%s line %d: %s", path, lineno, out);
            break;
        }
        if (out) defs++;
    }
    free(buf);
    fclose(f);
    memcpy(nm_on_error, saved, sizeof saved);
    if (err) nm_fail("%s", err);
    char *o = arena_alloc(len + 64);
    sprintf(o, "using %.*s (%d statements)", (int)len, name, defs);
    return o;
}

/* let f(x, y) = body   or   let A[0] = v0, A[1] = v1, A[k] = body */
static char *define_rule(const char *name, size_t len) {
    const char *def_start = toks[pos - 1].at;
    Def *d = perm_alloc(sizeof *d);
    memset(d, 0, sizeof *d);
    d->pn = perm_alloc(16 * sizeof(char *));
    Val v; memset(&v, 0, sizeof v);
    v.def = d;
    if (at_op('(')) {
        pos++;
        v.kind = V_FUNC;
        for (;;) {
            if (peek()->kind != T_NAME) nm_fail("expected a parameter name");
            if (d->np == 16) nm_fail("too many parameters");
            d->pn[d->np++] = dup_perm(peek()->s, peek()->len); pos++;
            if (at_op(',')) { pos++; continue; }
            expect_op(')');
            break;
        }
        expect_op('=');
        const char *bs = peek()->at;
        cases();
        if (peek()->kind != T_END) nm_fail("unexpected '%.*s' after the rule", (int)peek()->len, peek()->s);
        d->src = dup_perm(bs, (size_t)(peek()->at - bs));
    } else {
        v.kind = V_SEQ;
        d->iidx = perm_alloc(64 * sizeof(int64_t)); d->isrc = perm_alloc(64 * sizeof(char *));
        for (;;) {
            expect_op('[');
            if (peek()->kind == T_NUM) {                   /* a starting term */
                Val iv = vq(parse_number(peek()->s, peek()->len)); pos++;
                expect_op(']'); expect_op('=');
                if (d->ninit == 64) nm_fail("too many starting terms");
                d->iidx[d->ninit] = whole(iv, "an index");
                const char *bs = peek()->at;
                expr();
                d->isrc[d->ninit++] = dup_perm(bs, (size_t)(peek()->at - bs));
                expect_op(',');
                if (peek()->kind != T_NAME || !same(peek()->s, peek()->len, name, len))
                    nm_fail("expected the next term of %.*s", (int)len, name);
                pos++;
                continue;
            }
            if (peek()->kind != T_NAME) nm_fail("expected a starting term %.*s[0] = ... or the rule %.*s[k] = ...", (int)len, name, (int)len, name);
            d->pn[0] = dup_perm(peek()->s, peek()->len); d->np = 1; pos++;
            expect_op(']'); expect_op('=');
            const char *bs = peek()->at;
            cases();
            if (peek()->kind != T_END) nm_fail("unexpected '%.*s' after the rule", (int)peek()->len, peek()->s);
            d->src = dup_perm(bs, (size_t)(peek()->at - bs));
            break;
        }
    }
    bind(name, len, v);
    const char *def_end = peek()->at;
    while (def_end > def_start && (def_end[-1] == '\n' || def_end[-1] == '\r' || def_end[-1] == ' ')) def_end--;
    char *o = arena_alloc((size_t)(def_end - def_start) + 16);
    sprintf(o, "%.*s  [rule]", (int)(def_end - def_start), def_start);
    return o;
}

/* Run one statement; returns the line to print (arena memory), or NULL for an empty line. */
char *nm_run(const char *line, int *failed) {
    *failed = 0;
    if (setjmp(nm_on_error)) {
        *failed = 1;
        char *m = arena_alloc(strlen(nm_error_msg) + 8);
        sprintf(m, "error: %s", nm_error_msg);
        return m;
    }
    lex(line);
    pdepth = frame_base = calldepth = 0; nsum_skip = 0;
    stmt_id++;
    if (peek()->kind == T_END) return NULL;
    if (at_key("use")) {
        pos++;
        if (peek()->kind != T_NAME) nm_fail("expected a library name after 'use'");
        Tok t = *peek();
        char *name = dup_perm(t.s, t.len);             /* the library's own lines replace this one */
        return use_library(name, t.len);
    }
    const char *volatile let_name = NULL; volatile size_t let_len = 0;
    if (at_key("let")) {
        pos++;
        if (peek()->kind != T_NAME) nm_fail("expected a name after 'let'");
        let_name = peek()->s; let_len = peek()->len; pos++;
        if (at_op('(') || at_op('[')) return define_rule(let_name, let_len);
        expect_op('=');
    }
    const char *src_start = peek()->at;
    Node *volatile e = expr();
    const char *src_end = peek()->at;
    int64_t places = DEFAULT_PLACES, order = DEFAULT_ORDER;
    volatile int asked = 0;
    const char *volatile to_var = NULL; volatile size_t to_len = 0;
    if (at_key("to")) {
        pos++;
        if (peek()->kind == T_NAME) {                                  /* to x^8: x is the series letter */
            to_var = peek()->s; to_len = peek()->len;
            pos++; expect_op('^');
            if (peek()->kind != T_NUM || memchr(peek()->s, '.', peek()->len)) nm_fail("expected a whole number after '^'");
            Z k = z_from_dec(peek()->s, peek()->len);
            if (!z_fits_i64(k, &order) || order > 2000) nm_fail("order too large");
            pos++;
        } else {
            if (peek()->kind != T_NUM || memchr(peek()->s, '.', peek()->len)) nm_fail("expected a whole number after 'to'");
            Z k = z_from_dec(peek()->s, peek()->len);
            if (!z_fits_i64(k, &places) || places > 100000) nm_fail("too many places");
            pos++;
            expect_key("places");
            asked = 1;
        }
    }
    if (peek()->kind != T_END) nm_fail("unexpected '%.*s'", (int)peek()->len, peek()->s);
    work_prec = places + GUARD_DIGITS;
    char *text;
    Val v;
    if (e->kind == N_ROOT && e->para) {                /* Newton's parallelogram */
        if (let_name) nm_fail("a root found with the parallelogram cannot be named yet");
        return para_root(e, to_var, to_len, order);
    }
    if (!setjmp(need_series)) {
        v = eval(e);
        text = show(v, places, asked);
        if (v.kind == V_BALL) {                        /* keep the recipe, so a later request can carry it further */
            memset(&v, 0, sizeof v);
            v.kind = V_NUMREC;
            v.src = dup_perm(src_start, (size_t)(src_end - src_start));
        }
    } else {
        const char *var = to_var; size_t vlen = to_len;
        if (!var) {
            Letters L; L.n = 0;
            letters(e, &L, NULL, 0);
            int nf = 0, xf = -1;
            for (int i = 0; i < L.n; i++) if (flowing_letter(L.s[i], L.len[i])) { nf++; xf = i; }
            if (L.n > 1 && nf == 1) { var = L.s[xf]; vlen = L.len[xf]; }   /* Newton: a, b, c are given; x, z flow */
            else if (L.n > 1) nm_fail("several letters (%.*s, %.*s ...): name the series letter with 'to %.*s^8'",
                                      (int)L.len[0], L.s[0], (int)L.len[1], L.s[1], (int)L.len[1], L.s[1]);
            else if (L.n == 1) { var = L.s[0]; vlen = L.len[0]; }
            else { var = "x"; vlen = 1; }
        }
        char *vname = dup_perm(var, vlen);
        SCtx cx; memset(&cx, 0, sizeof cx);
        cx.var = vname; cx.n = (int)order + 1;
        Dual d = sev(e, &cx);
        text = s_to_str(d.a, vname, (int)order);
        memset(&v, 0, sizeof v);
        v.kind = V_REC;
        v.src = dup_perm(src_start, (size_t)(src_end - src_start));
        v.var = vname;
    }
    if (let_name && v.kind == V_ROOT) {                   /* a named root: a letter with its equation */
        C r = c_surd_from_root(dup_perm(let_name, let_len), v.root);
        Q k;
        if (c_const_value(r, &k)) v = vq(k); else { memset(&v, 0, sizeof v); v.kind = V_POLY; v.c = r; }
    }
    if (let_name) {
        bind(let_name, let_len, v);
        char *o = arena_alloc(strlen(text) + let_len + 4);
        sprintf(o, "%.*s = %s", (int)let_len, let_name, text);
        return o;
    }
    return text;
}
