/* The language: reading a statement, working it out, and writing the result.
 *
 *   statement := "use" NAME | "let" NAME "=" expr [to] | expr [to]
 *   to        := "to" INTEGER "places" | "to" LETTER "^" INTEGER
 *   expr      := term { ("+" | "-") term }
 *   term      := unary { ("*" | "/") unary | unary }          (side by side multiplies: 2y, 3(x+1))
 *   unary     := "-" unary | power
 *   power     := primary [ "^" unary ]
 *   primary   := NUMBER | NAME{'} | NAME "(" expr ")" | "(" expr ")" | "sqrt" "(" expr ")"
 *              | "d/d"LETTER primary | "integral" "(" expr ["," LETTER ["," expr "," expr]] ")"
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

static const char *KEYWORDS[] = {"let", "to", "places", "root", "of", "near", "sqrt", "integral", "use", "for", "starting", "sum", "if", "otherwise", "solve", "eliminate", "from", "show", NULL};

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

static int same_word(Tok t, const char *w) { return t.len == strlen(w) && !strncmp(t.s, w, t.len); }
static int reads_as_letters(const char *s, size_t len);
static Ball conic_pi(void *data, int64_t prec);
static int number_sign(C a);
static C conic_constant(int kind, C arg);
static void root_legend(C *vals, int nv, char **o);

static void lex(const char *src) {
    size_t n = strlen(src);
    toks = arena_alloc((n * 3 + 2) * sizeof(Tok));
    ntok = 0; pos = 0;
    int header_open = 0, nheader = 0;
    const char *header_s[16]; size_t header_l[16];
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
            int key = is_key(src + i, j - i);
            size_t k = j;
            while (k < n && isspace((unsigned char)src[k])) k++;
            int exempt = key || j - i == 1 || (k < n && strchr("(['", src[k]))
                         || (ntok >= 1 && toks[ntok - 1].kind == T_KEY && (same_word(toks[ntok - 1], "let") || same_word(toks[ntok - 1], "use")))
                         || (ntok >= 2 && toks[ntok - 1].kind == T_OP && toks[ntok - 1].op == '/' && same_word(toks[ntok - 2], "d"))
                         || header_open;
            if (header_open && !key && nheader < 16) { header_s[nheader] = src + i; header_l[nheader++] = j - i; }
            for (int h = 0; h < nheader && !exempt; h++) exempt = header_l[h] == j - i && !strncmp(header_s[h], src + i, j - i);
            if (!exempt && reads_as_letters(src + i, j - i)) {          /* Newton's ax: the product of a and x */
                for (size_t q = i; q < j; q++) push(T_NAME, src + q, 1, 0);
            } else push(key ? T_KEY : T_NAME, src + i, j - i, 0);
            while (j < n && src[j] == '\'') { toks[ntok - 1].primes++; j++; }   /* y', y'' */
            if (ntok == 2 && same_word(toks[0], "let") && j < n && src[j] == '(') header_open = 1;   /* let f(a, b): names */
            i = j; continue;
        }
        if (header_open && c == ')') {
            header_open = 0;
            push(T_OP, src + i, 1, ')'); i++; continue;
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

enum { N_NUM, N_NAME, N_NEG, N_BIN, N_SQRT, N_ROOT, N_APPLY, N_DERIV, N_INTEG, N_INDEX, N_SUM, N_CASES, N_CMP, N_MAT };

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

static int name_known(const char *s, size_t len);
static Node *expr(void);
static Node *unary(void);
static Node *power(void);

static Node *primary(void) {
    Tok *t = peek();
    if (t->kind == T_NUM) { Node *n = mk(N_NUM); n->s = t->s; n->len = t->len; pos++; return n; }
    if (t->kind == T_NAME) {
        Tok *t1 = peek_at(1), *t2 = peek_at(2);
        if (t->len == 1 && !t->primes && t1->kind == T_OP && t1->op == '(' && !(t->s[0] == 'd' && t2->kind == T_OP)
            && !name_known(t->s, 1)) {          /* x(x + 1)^2: a letter times, as Newton writes */
            Node *n = mk(N_NAME); n->s = t->s; n->len = 1; pos++;
            return n;
        }
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
    if (at_op('[')) {                                                  /* [[1, 2], [3, 4]], or [1, 2, 3] (a column) */
        pos++;
        Node *n = mk(N_MAT);
        n->args = arena_alloc(256 * sizeof(Node *));
        if (at_op('[')) {
            while (at_op('[')) {
                pos++;
                Node *row = mk(N_MAT);
                row->args = arena_alloc(256 * sizeof(Node *));
                do { if (row->nargs == 256) nm_fail("too many entries"); row->args[row->nargs++] = expr(); } while (at_op(',') && (pos++, 1));
                expect_op(']');
                if (n->nargs == 256) nm_fail("too many rows");
                n->args[n->nargs++] = row;
                if (at_op(',')) pos++;
            }
            n->op = 'm';
        } else {
            do { if (n->nargs == 256) nm_fail("too many entries"); n->args[n->nargs++] = expr(); } while (at_op(',') && (pos++, 1));
            n->op = 'v';
        }
        expect_op(']');
        return n;
    }
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
            if (at_op(',')) { pos++; n->b=expr(); expect_op(','); n->c=expr(); }
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

enum { V_Q, V_POLY, V_ROOT, V_BALL, V_REC, V_NUMREC, V_FUNC, V_SEQ, V_MAT, V_CBALL, V_TEXT, V_RAT, V_AREA, V_APART };
struct Def;
typedef struct {
    int kind; Q q; C c; Root *root; Ball ball; const char *src, *var; struct Def *def;
    R rat; Integral area; Apart apart; const char *verdict;
    Mat m; CBall z; const char *text;       /* a matrix; a complex approximate value; a verdict to print */
} Val;   /* V_POLY holds c */

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
static Val vr(R r) {
    Q k;
    if (c_const_value(r.den, &k)) return vc(c_scale(r.num, q_div(q_from_z(z_from_i64(1)), k)));
    Val v; memset(&v, 0, sizeof v); v.kind = V_RAT; v.rat = r; return v;
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

static C as_c(Val v) {
    if (v.kind == V_APART) return as_c(vr(apart_sum(v.apart)));
    if (v.kind == V_POLY) return v.c;
    if (v.kind == V_Q) return c_const(v.q);
    if (v.kind == V_RAT && c_is_monomial(v.rat.den)) return c_div(v.rat.num, v.rat.den);
    nm_fail("a polynomial quantity is required here; a rational function has a denominator");
}
static R as_r(Val v) { return v.kind == V_APART ? apart_sum(v.apart) : v.kind == V_RAT ? v.rat : r_from_c(as_c(v)); }

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

static int is_exact(Val v) { return v.kind == V_Q || v.kind == V_POLY || v.kind == V_RAT || v.kind == V_APART; }

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
    if (letter_is_named(l)) return named_ball(l,prec);
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
            if (!q_sign(ct_e(&a.t[t], l))) continue;
            if (!z_fits_i64(ct_e(&a.t[t], l).num, &e) || !q_is_int(ct_e(&a.t[t], l))) nm_fail("internal: a fractional power of a surd");
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
        if (li >= 0 && (!q_is_int(ct_e(&c.t[t], li)) || !z_fits_i64(ct_e(&c.t[t], li).num, &e) || e < 0 || e > 100000))
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
        if (li >= 0) z_fits_i64(ct_e(&c.t[t], li).num, &e);
        p.c[e] = q_add(p.c[e], c.t[t].k);
    }
    while (p.deg >= 0 && q_sign(p.c[p.deg]) == 0) p.deg--;
    return p;
}

/* ---------------- the exact pass: numbers, rationals, polynomials ---------------- */

static Val eval(Node *n);

/* while reading an equation's coefficients, the unknown and its derivatives stand for given numbers */
static struct { int active; const char *u; size_t ul; Q vals[8]; } ovr;

/* inside integral(..): the root of the flowing quantity is the letter put for it (sqrtint.c) */
static struct { int on; C q; int s; } sqo;

static Val power_val(Val base, Val ex) {
    if (sqo.on && ex.kind == V_Q && !z_is_one(ex.q.den) && z_cmp(ex.q.den, z_from_i64(2)) == 0
        && (base.kind == V_Q || base.kind == V_POLY) && c_equal(as_c(base), sqo.q)) {
        int64_t k; z_fits_i64(ex.q.num, &k);
        return vr(r_pow_int(r_from_c(c_letter(sqo.s)), k));
    }
    if (ex.kind != V_Q) nm_fail("the exponent must be an exact number");
    Q e = ex.q;
    int64_t num, den;
    if (!z_fits_i64(e.num, &num) || !z_fits_i64(e.den, &den) || num > 1000000 || num < -1000000)
        nm_fail("exponent too large");
    if (base.kind == V_Q) {
        if (den == 1) return vq(q_pow(base.q, num));
        return nth_root(q_pow(base.q, num), den);
    }
    if (base.kind == V_RAT) {
        if (den != 1) nm_fail("fractional powers of rational functions are not supported");
        return vr(r_pow_int(base.rat, num));
    }
    if (base.kind == V_POLY) {
        if (den == 1 && num >= 0) return vc(c_pow_int(base.c, num));
        if (den == 1 && num < 0) return vr(r_pow_int(r_from_c(base.c), num));
        C out;
        if (den != 1 && c_pow_q(base.c, e, &out)) return vc(out);
        longjmp(need_series, 1);
    }
    if (den != 1) nm_fail("fractional powers of approximate numbers are not in this version");
    return vball(b_pow(as_ball(base), num, work_prec));
}

static int rational_letter(R r) {
    int v = -1;
    for (int l = 0; l < letter_count(); l++) if (c_uses(r.num,l) || c_uses(r.den,l)) {
        if (letter_is_surd(l) || letter_is_named(l)) continue;
        if (v >= 0) return -1;
        v = l;
    }
    return v;
}
static R substitute_r(C p, int v, R x) {
    UP a = up_from(p,v);
    R out = r_from_c(c_zero());
    for (int i = a.deg; i >= 0; i--) out = r_add(r_mul(out,x),r_from_c(a.c[i]));
    return out;
}
static C derivative_c(C p, int l) {
    C out = c_zero();
    for (int t = 0; t < p.nt; t++) {
        CT term = p.t[t];
        if (!q_sign(ct_e(&term, l))) continue;
        term.k = q_mul(term.k,ct_e(&term, l)); ct_set(&term, l, q_sub(ct_e(&term, l),qi(1)));
        C m = {1,&term}; out = c_add(out,m);
    }
    return out;
}

static int max_primes(Node *n, const char *u, size_t ul);
static Val vmat(Mat m) { Val v; memset(&v, 0, sizeof v); v.kind = V_MAT; v.m = m; return v; }
static Val vcb(CBall z) { Val v; memset(&v, 0, sizeof v); v.kind = V_CBALL; v.z = z; return v; }
static Val vtext(const char *t) { Val v; memset(&v, 0, sizeof v); v.kind = V_TEXT; v.text = t; return v; }
static int split_imag(C a, C *re, C *im);
static Ball ball_of_c(C a, int64_t prec);

static CBall as_cball(Val v) {
    CBall z;
    if (v.kind == V_CBALL) return v.z;
    if (v.kind == V_POLY) {
        C re, im;
        if (split_imag(v.c, &re, &im)) { z.re = ball_of_c(re, work_prec); z.im = ball_of_c(im, work_prec); return z; }
    }
    if (v.kind == V_MAT || v.kind == V_TEXT) nm_fail("this is not a number");
    z.re = as_ball(v); z.im = b_from_q(qi(0), work_prec);
    return z;
}

static R entry_of(Val v) {                              /* a matrix entry, or a scalar for a matrix */
    if (is_exact(v)) return as_r(v);
    nm_fail("matrix entries must be exact (numbers, surds or letters)");
}

static Val arith_value(char op, Val a, Val b) {
    if (a.kind==V_AREA || b.kind==V_AREA) nm_fail("general arithmetic on conic areas comes later; differentiate or substitute a number");
    if (a.kind==V_APART) a=vr(apart_sum(a.apart));
    if (b.kind==V_APART) b=vr(apart_sum(b.apart));
    if (a.kind == V_TEXT || b.kind == V_TEXT) nm_fail("a verdict cannot enter a calculation");
    if (a.kind == V_MAT || b.kind == V_MAT) {
        if (op == '^') {
            if (a.kind != V_MAT || b.kind != V_Q || !q_is_int(b.q)) nm_fail("a matrix power needs a whole exponent");
            int64_t e; z_fits_i64(b.q.num, &e);
            Mat m = mat_pow(a.m,e);
            Val v = vmat(m);
            if (e < 0) v.verdict = mat_verdict(a.m,m,-1);
            return v;
        }
        if (a.kind == V_MAT && b.kind == V_MAT) {
            if (op == '+') return vmat(mat_add(a.m, b.m, 1));
            if (op == '-') return vmat(mat_add(a.m, b.m, -1));
            if (op == '*') return vmat(mat_mul(a.m, b.m));
            nm_fail("dividing by a matrix: multiply by its inverse, M^-1");
        }
        if (op == '*') return a.kind == V_MAT ? vmat(mat_scale(a.m, entry_of(b))) : vmat(mat_scale(b.m, entry_of(a)));
        if (op == '/' && a.kind == V_MAT) return vmat(mat_scale(a.m, r_div(r_from_c(c_const(qi(1))), entry_of(b))));
        nm_fail("a matrix and a number can only be multiplied, or the matrix divided by the number");
    }
    if (a.kind == V_CBALL || b.kind == V_CBALL) {
        if (op == '^') {
            if (b.kind != V_Q || !q_is_int(b.q)) nm_fail("a complex approximate value takes whole powers only");
            int64_t e; z_fits_i64(b.q.num, &e);
            CBall r = cb_from_q(qi(1), qi(0), work_prec), x = as_cball(a);
            int neg = e < 0; if (neg) e = -e;
            while (e) { if (e & 1) r = cb_mul(r, x, work_prec); e >>= 1; if (e) x = cb_mul(x, x, work_prec); }
            return vcb(neg ? cb_div(cb_from_q(qi(1), qi(0), work_prec), r, work_prec) : r);
        }
        CBall x = as_cball(a), y = as_cball(b);
        switch (op) {
        case '+': return vcb(cb_add(x, y, work_prec));
        case '-': return vcb(cb_sub(x, y, work_prec));
        case '*': return vcb(cb_mul(x, y, work_prec));
        case '/': return vcb(cb_div(x, y, work_prec));
        }
    }
    if (op == '^') return power_val(a, b);
    if (a.kind == V_Q && b.kind == V_Q) {
        switch (op) {
        case '+': return vq(q_add(a.q, b.q));
        case '-': return vq(q_sub(a.q, b.q));
        case '*': return vq(q_mul(a.q, b.q));
        case '/': return vq(q_div(a.q, b.q));
        }
    }
    if (op == '/' && is_exact(a) && a.kind != V_RAT && (b.kind == V_Q || (b.kind == V_POLY && !c_has_plain(b.c)))) {
        C d = as_c(b);                                /* by a number or surd: no quotient of letters is needed */
        if (c_is_zero(d)) nm_fail("division by zero");
        return vc(c_mul(as_c(a), c_inv(d)));
    }
    if (is_exact(a) && is_exact(b) && (a.kind == V_RAT || b.kind == V_RAT || op == '/')) {
        R x = as_r(a), y = as_r(b);
        switch (op) {
        case '+': return vr(r_add(x, y));
        case '-': return vr(r_sub(x, y));
        case '*': return vr(r_mul(x, y));
        case '/': return vr(r_div(x, y));
        }
    }
    if (is_exact(a) && is_exact(b)) {
        C x = as_c(a), y = as_c(b);
        switch (op) {
        case '+': return vc(c_add(x, y));
        case '-': return vc(c_sub(x, y));
        case '*': return vc(c_mul(x, y));
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

static const char *combine_verdicts(const char *a, const char *b) {
    if (a && !strcmp(a,"[exact]")) a = NULL;
    if (b && !strcmp(b,"[exact]")) b = NULL;
    if (!a) return b;
    if (!b || !strcmp(a,b)) return a;
    char *out = arena_alloc(strlen(a) + strlen(b) + 48);
    sprintf(out,"[exact; subject to %s and %s]",a,b);
    return out;
}
static Val arith(char op, Val a, Val b) {
    Val v = arith_value(op,a,b);
    v.verdict = combine_verdicts(v.verdict,combine_verdicts(a.verdict,b.verdict));
    return v;
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
    if (same(s, len, "pi", 2)) return vc(c_named("pi", conic_pi, NULL));
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

static Val series_value_d(Binding *b, Val arg, int deriv);
static C conic_constant(int kind,C arg);
static void integral_no_poles(R f,int v,C lo,C hi);
static int recipe_is_equation(Binding *b);

static Z whole_z(Val v, const char *what) {
    if (v.kind != V_Q || !q_is_int(v.q)) nm_fail("%s must be a whole number", what);
    return v.q.num;
}

static Mat as_mat(Val v, const char *what) {
    if (v.kind != V_MAT) nm_fail("%s must be a matrix", what);
    return v.m;
}

static char *verdict_of_prime(int s) { return s == 1 ? "[proved]" : s == 2 ? "[probable: BPSW test, no proof found]" : "[exact]"; }

static int inherit_builtin(Val *out, Val *args, int n) {
    const char *conditions = NULL;
    for (int i = 0; i < n; i++) conditions = combine_verdicts(conditions,args[i].verdict);
    if (conditions && out->kind == V_TEXT) {
        char *s = arena_alloc(strlen(out->text) + strlen(conditions) + 24);
        sprintf(s,"%s [subject to %s]",out->text,conditions); out->text = s;
    } else out->verdict = combine_verdicts(out->verdict,conditions);
    return 1;
}

/* the functions every mathematician expects, when the name is not defined otherwise */
static const char *BUILTINS[] = {"det", "inverse", "transpose", "rank", "nullspace", "charpoly", "eigenvalues", "linsolve",
    "gcd", "lcm", "mod", "powmod", "invmod", "isprime", "apart", "factor", "divisors", "sigma", "phi", "nextprime",
    "log", "atan", "pi", NULL};

static int name_known(const char *s, size_t len) { return lookup(s, len) || param_find(s, len); }

/* an unbound name of letters only, not a function, parameter or defined name, is a product of single letters */
static int reads_as_letters(const char *s, size_t len) {
    for (size_t i = 0; i < len; i++) if (!isalpha((unsigned char)s[i])) return 0;
    for (int i = 0; BUILTINS[i]; i++) if (same(s, len, BUILTINS[i], strlen(BUILTINS[i]))) return 0;
    return !lookup(s, len) && !param_find(s, len);
}

static int builtin(const char *s, size_t len, Node *n, Val *out) {
#define IS(name) same(s, len, name, strlen(name))
    Val a[8];
    int na = n->nargs;
    if (na > 8) return 0;
    if (!(IS("det") || IS("inverse") || IS("transpose") || IS("rank") || IS("nullspace") || IS("charpoly") || IS("eigenvalues")
          || IS("linsolve") || IS("gcd") || IS("lcm") || IS("mod") || IS("powmod") || IS("invmod") || IS("isprime")
          || IS("apart") || IS("factor") || IS("divisors") || IS("sigma") || IS("phi") || IS("nextprime")
          || IS("log") || IS("atan"))) return 0;
    for (int i = 0; i < na; i++) a[i] = eval(n->args[i]);
    if (IS("log") || IS("atan")) {                    /* areas of the hyperbola and the circle, as exact constants */
        if (na != 1) nm_fail("%.*s takes one number", (int)len, s);
        C u;
        if (a[0].kind == V_Q) u = c_const(a[0].q);
        else if (a[0].kind == V_POLY && !c_has_plain(a[0].c)) u = a[0].c;
        else nm_fail("%.*s here takes an exact number; for a series use log1p or atan from the prelude", (int)len, s);
        if (IS("log") && number_sign(u) <= 0) nm_fail("log of a number that is not positive");
        C r = conic_constant(IS("log") ? AREA_LOG : AREA_ATAN, u);
        Q k;
        *out = c_const_value(r, &k) ? vq(k) : vc(r);
        return 1;
    }
#define NEED(k) do { if (na != (k)) nm_fail("%.*s takes %d argument%s", (int)len, s, (k), (k) == 1 ? "" : "s"); } while (0)
    if (IS("apart")) {
        NEED(2);
        if(n->args[1]->kind!=N_NAME || n->args[1]->primes) nm_fail("apart needs its integration letter as the second argument");
        *out=vq(qi(0));out->kind=V_APART;
        out->apart=integ_apart(as_r(a[0]),letter_index(n->args[1]->s,n->args[1]->len));
        return inherit_builtin(out,a,na);
    }
    if (IS("det")) { NEED(1); *out = vr(mat_det(as_mat(a[0], "the argument of det"))); return inherit_builtin(out,a,na); }
    if (IS("inverse")) {
        NEED(1); Mat m = as_mat(a[0], "the argument of inverse"), inv = mat_inverse(m);
        *out = vmat(inv); out->verdict = mat_verdict(m, inv, -1); return inherit_builtin(out,a,na);
    }
    if (IS("transpose")) { NEED(1); *out = vmat(mat_transpose(as_mat(a[0], "the argument of transpose"))); return inherit_builtin(out,a,na); }
    if (IS("rank")) {
        NEED(1); Mat m = as_mat(a[0], "the argument of rank"); int rank = mat_rank(m);
        *out = vq(qi(rank)); out->verdict = mat_verdict(m, mat_new(0,0), rank); return inherit_builtin(out,a,na);
    }
    if (IS("nullspace")) {
        NEED(1); Mat m = as_mat(a[0], "the argument of nullspace"), ns = mat_nullspace(m);
        *out = vmat(ns); out->verdict = mat_verdict(m, ns, m.c - ns.c); return inherit_builtin(out,a,na);
    }
    if (IS("charpoly") || IS("eigenvalues")) {
        if (na != 1 && !(IS("charpoly") && na == 2)) nm_fail("%.*s takes a matrix", (int)len, s);
        Mat m = as_mat(a[0], "the argument");
        int t;
        if (na == 2) {
            if (n->args[1]->kind != N_NAME) nm_fail("the second argument of charpoly is a letter");
            t = letter_index(n->args[1]->s, n->args[1]->len);
        } else t = letter_index("t", 1);
        R *cp = mat_charpoly(m), p = r_from_c(c_zero());
        for (int i = 0; i <= m.r; i++) p = r_add(p, r_mul(cp[i], r_from_c(c_pow_int(c_letter(t), i))));
        if (IS("charpoly")) { *out = vr(p); return inherit_builtin(out,a,na); }
        C roots[256]; int unres;
        if (c_uses(p.den, t)) nm_fail("the characteristic letter occurs in an entry denominator");
        int k = elim_roots(p.num, t, roots, 256, &unres);
        Mat r = mat_new(1, k);
        for (int i = 0; i < k; i++) r.a[i] = r_from_c(roots[i]);
        if (!unres) { *out = vmat(r); return inherit_builtin(out,a,na); }
        char *txt = arena_alloc(strlen(mat_to_str(r)) + 128);
        sprintf(txt, "%s  [and %d complex roots of a factor above degree two, not shown]", mat_to_str(r), unres);
        *out = vtext(txt);
        return inherit_builtin(out,a,na);
    }
    if (IS("linsolve")) {
        NEED(2);
        Mat ns, x = mat_solve(as_mat(a[0], "the matrix"), as_mat(a[1], "the right side"), &ns);
        Mat m = as_mat(a[0], "the matrix");
        Mat formulas = mat_new(x.r, x.c + ns.c);
        for (int i = 0; i < x.r; i++) {
            for (int j = 0; j < x.c; j++) formulas.a[i * formulas.c + j] = x.a[i * x.c + j];
            for (int j = 0; j < ns.c; j++) formulas.a[i * formulas.c + x.c + j] = ns.a[i * ns.c + j];
        }
        const char *verdict = mat_verdict(m, formulas, m.c - ns.c);
        if (ns.c == 0) { *out = vmat(x); out->verdict = verdict; return inherit_builtin(out,a,na); }
        char *t1 = mat_to_str(x), *t2 = mat_to_str(ns);
        char *txt = arena_alloc(strlen(t1) + strlen(t2) + strlen(verdict) + 96);
        sprintf(txt, "%s + any combination of the columns of %s  %s", t1, t2, verdict);
        *out = vtext(txt);
        return inherit_builtin(out,a,na);
    }
    if (IS("gcd")) {
        int polynomials = 0;
        for (int i = 0; i < na; i++) if (a[i].kind == V_POLY) polynomials = 1;
        if (polynomials) {
            C g = c_zero();
            for (int i = 0; i < na; i++) g = poly_gcd(g, as_c(a[i]));
            *out = vc(g); return inherit_builtin(out,a,na);
        }
    }
    if (IS("gcd") || IS("lcm")) {
        if (na < 1) nm_fail("%.*s needs numbers", (int)len, s);
        Z g = whole_z(a[0], "an argument");
        g = z_abs(g);
        for (int i = 1; i < na; i++) {
            Z b = z_abs(whole_z(a[i], "an argument"));
            if (IS("gcd")) g = z_gcd(g, b);
            else { Z q, r, d = z_gcd(g, b); if (d.s) { z_divmod(z_mul(g, b), d, &q, &r); g = q; } else g = z_zero(); }
        }
        *out = vq(q_from_z(g));
        return inherit_builtin(out,a,na);
    }
    if (IS("mod")) { NEED(2); *out = vq(q_from_z(z_mod(whole_z(a[0], "a"), whole_z(a[1], "m")))); return inherit_builtin(out,a,na); }
    if (IS("powmod")) {
        NEED(3);
        Z b = whole_z(a[0], "the base"), e = whole_z(a[1], "the exponent"), m = whole_z(a[2], "the modulus");
        if (e.s < 0) { Z iv; if (!z_invmod(b, m, &iv)) nm_fail("%s has no inverse modulo %s", z_to_str(b), z_to_str(m)); b = iv; e = z_neg(e); }
        *out = vq(q_from_z(z_powmod(b, e, m)));
        return inherit_builtin(out,a,na);
    }
    if (IS("invmod")) {
        NEED(2);
        Z iv, b = whole_z(a[0], "a"), m = whole_z(a[1], "m");
        if (!z_invmod(b, m, &iv)) nm_fail("%s has no inverse modulo %s (they share a factor)", z_to_str(b), z_to_str(m));
        *out = vq(q_from_z(iv));
        return inherit_builtin(out,a,na);
    }
    if (IS("isprime")) {
        NEED(1);
        int v = z_isprime(whole_z(a[0], "the argument"));
        char *txt = arena_alloc(96);
        sprintf(txt, "%s  %s", v ? "true" : "false", verdict_of_prime(v));
        *out = vtext(txt);
        return inherit_builtin(out,a,na);
    }
    if (IS("factor") && na == 1 && a[0].kind == V_POLY) {     /* a polynomial: Newton's divisors */
        int nl = 0;
        for (int l = 0; l < letter_count(); l++) if (c_uses(a[0].c, l)) nl++;
        if (nl != 1) nm_fail("factor: a polynomial in one letter (several letters: later)");
        PolyFactors P = poly_factor(c_to_poly(a[0].c, NULL));
        size_t cap = 128;
        for (int i = 0; i < P.n; i++) cap += strlen(p_to_str(P.f[i])) + 32;
        cap += strlen(q_to_str(P.unit));
        char *txt = arena_alloc(cap), *o = txt;
        if (q_cmp(P.unit, qi(-1)) == 0) o += sprintf(o, "-");
        else if (q_cmp(P.unit, qi(1))) o += sprintf(o, q_is_int(P.unit) ? "%s" : "(%s)", q_to_str(P.unit));
        for (int i = 0; i < P.n; i++) {
            int terms = 0;
            for (int j = 0; j <= P.f[i].deg; j++) terms += q_sign(P.f[i].c[j]) != 0;
            o += sprintf(o, terms == 1 || (P.n == 1 && P.e[i] == 1 && !q_cmp_one(P.unit)) ? "%s" : "(%s)", p_to_str(P.f[i]));
            if (P.e[i] > 1) o += sprintf(o, "^%d", P.e[i]);
        }
        sprintf(o, "  [multiplied back; every factor proved irreducible over Q (all recombinations tried)]");
        {   /* facts read off the factors: nothing more is computed */
            size_t fc = 64;
            for (int i = 0; i < P.n; i++) fc += 48 + (P.f[i].deg == 1 ? strlen(q_to_str(P.f[i].c[0])) + strlen(q_to_str(P.f[i].c[1])) + 8 : 0);
            char *deg = arena_alloc(fc), *mul = arena_alloc(fc), *roots = arena_alloc(fc), *d = deg, *m = mul, *r = roots;
            int total = 0, sqfree = 1;
            d += sprintf(d, "["); m += sprintf(m, "["); r += sprintf(r, "[");
            for (int i = 0; i < P.n; i++) {
                d += sprintf(d, "%s%d", i ? "," : "", P.f[i].deg);
                m += sprintf(m, "%s%d", i ? "," : "", P.e[i]);
                total += P.f[i].deg * P.e[i];
                if (P.e[i] > 1) sqfree = 0;
                if (P.f[i].deg == 1) r += sprintf(r, "%s%s", r[-1] == '[' ? "" : ",", nm_json_str(q_to_str(q_neg(q_div(P.f[i].c[0], P.f[i].c[1])))));
            }
            sprintf(d, "]"); sprintf(m, "]"); sprintf(r, "]");
            nm_fact("degree", "%d", total);
            nm_fact("irreducible", P.n == 1 && P.e[0] == 1 ? "true" : "false");
            nm_fact("squarefree", sqfree ? "true" : "false");
            nm_fact("factor_degrees", "%s", deg);
            nm_fact("multiplicities", "%s", mul);
            nm_fact("rational_roots", "%s", roots);
        }
        *out = vtext(txt);
        return inherit_builtin(out,a,na);
    }
    if (IS("factor")) {
        NEED(1);
        Factors F = z_factor(whole_z(a[0], "the argument"));
        size_t cap = 128;
        for (int i = 0; i < F.n; i++) cap += strlen(z_to_str(F.p[i])) + 32;
        char *txt = arena_alloc(cap), *o = txt;
        int probable = 0;
        if (F.sign < 0) o += sprintf(o, "-");
        if (F.n == 0) o += sprintf(o, "1");
        for (int i = 0; i < F.n; i++) {
            o += sprintf(o, "%s%s", i ? " * " : "", z_to_str(F.p[i]));
            if (F.e[i] > 1) o += sprintf(o, "^%d", F.e[i]);
            if (F.status[i] == 2) probable++;
        }
        if (F.n == 0) sprintf(o, "  [exact]");
        else if (probable) sprintf(o, "  [multiplied back; %d factor%s only probable prime%s (BPSW)]", probable, probable > 1 ? "s" : "", probable > 1 ? "s" : "");
        else sprintf(o, "  [multiplied back; every factor proved prime]");
        nm_fact("prime", F.n == 1 && F.e[0] == 1 && F.sign > 0 ? (F.status[0] == 1 ? "true" : "\"probable\"") : "false");
        nm_fact("distinct_prime_factors", "%d", F.n);
        *out = vtext(txt);
        return inherit_builtin(out,a,na);
    }
    if (IS("divisors")) {
        NEED(1);
        Z *dv = arena_alloc(100000 * sizeof(Z));
        int k = z_divisors(z_abs(whole_z(a[0], "the argument")), dv, 100000);
        Mat r = mat_new(1, k);
        for (int i = 0; i < k; i++) r.a[i] = r_from_c(c_const(q_from_z(dv[i])));
        *out = vmat(r);
        return inherit_builtin(out,a,na);
    }
    if (IS("sigma")) { NEED(1); *out = vq(q_from_z(z_sigma(z_abs(whole_z(a[0], "the argument"))))); return inherit_builtin(out,a,na); }
    if (IS("phi")) { NEED(1); *out = vq(q_from_z(z_phi(z_abs(whole_z(a[0], "the argument"))))); return inherit_builtin(out,a,na); }
    if (IS("nextprime")) { NEED(1); *out = vq(q_from_z(z_nextprime(whole_z(a[0], "the argument")))); return inherit_builtin(out,a,na); }
#undef IS
#undef NEED
    return 0;
}

static int mentions(Node *n, const char *v, size_t vl) {
    if (!n) return 0;
    if (n->kind == N_NAME && same(n->s, n->len, v, vl)) return 1;
    if (mentions(n->a, v, vl) || mentions(n->b, v, vl) || mentions(n->c, v, vl)) return 1;
    for (int i = 0; i < n->nargs; i++) if (mentions(n->args[i], v, vl)) return 1;
    return 0;
}

/* the radicands under sqrt(..) or ^(k/2) that hold the flowing letter; 0 none, 1 one, 2 several */
static int find_root(Node *n, const char *v, size_t vl, Node **rad) {
    if (!n) return 0;
    Node *cand = NULL;
    if (n->kind == N_SQRT && mentions(n->a, v, vl)) cand = n->a;
    if (n->kind == N_BIN && n->op == '^' && mentions(n->a, v, vl) && !mentions(n->b, v, vl)) {
        Val e = eval(n->b);
        if (e.kind == V_Q && z_cmp(e.q.den, z_from_i64(2)) == 0) cand = n->a;
    }
    int r = 0;
    if (cand) {
        if (*rad && !c_equal(as_c(eval(*rad)), as_c(eval(cand)))) return 2;
        *rad = cand; r = 1;
    }
    Node *kids[3] = {n->a, n->b, n->c};
    for (int i = 0; i < 3; i++) { int k = find_root(kids[i], v, vl, rad); if (k == 2) return 2; if (k) r = 1; }
    for (int i = 0; i < n->nargs; i++) { int k = find_root(n->args[i], v, vl, rad); if (k == 2) return 2; if (k) r = 1; }
    return r;
}

static C conic_constant(int kind, C arg);

/* integral(f(x, sqrt(q)), x [, lo, hi]): Newton's letter for the root, then the conic areas */
static int root_integral(Node *n, Val *out) {
    const char *v = n->len ? n->s : "x"; size_t vl = n->len ? n->len : 1;
    Node *rad = NULL;
    int k = find_root(n->a, v, vl, &rad);
    if (!k) return 0;
    if (k == 2) nm_fail("several different roots in one integrand come later");
    Val qv = eval(rad);
    if (qv.kind != V_POLY && qv.kind != V_Q) nm_fail("the quantity under the root must be a polynomial in %.*s", (int)vl, v);
    C q = as_c(qv);
    if (q.nt == 1) return 0;                          /* x^(3/2): a power of a letter is already exact */
    int x = letter_index(v, vl);
    char *disp = arena_alloc(strlen(c_to_str(q)) + 8);
    sprintf(disp, "sqrt(%s)", c_to_str(q));
    int s = letter_shown(disp);
    sqo.on = 1; sqo.q = q; sqo.s = s;
    jmp_buf saved; memcpy(saved, need_series, sizeof saved);
    if (setjmp(need_series)) { sqo.on = 0; memcpy(need_series, saved, sizeof saved); nm_fail("this integrand has no finite area here; ask for the series with 'to %.*s^N'", (int)vl, v); }
    Val f = eval(n->a);
    memcpy(need_series, saved, sizeof saved);
    sqo.on = 0;
    if (!is_exact(f)) nm_fail("the integrand must be exact");
    R fr = as_r(f);
    Integral I = sqrt_integral(fr, x, s, q);
    if (n->b) {
        Val lv = eval(n->b), hv = eval(n->c);
        if (!is_exact(lv) || !is_exact(hv)) nm_fail("definite integration needs exact real endpoints");
        C lo = as_c(lv), hi = as_c(hv);
        R A, B;                                       /* poles: all of A's; B's except a simple root of q, */
        sqrt_integral_split(fr, s, q, &A, &B);        /* where B s ~ 1/sqrt(q) keeps a finite area */
        if (c_uses(A.den, x)) integral_no_poles(r_make(c_const(qi(1)), A.den), x, lo, hi);
        R Bq = r_mul(B, r_from_c(q));
        if (c_uses(Bq.den, x)) integral_no_poles(r_make(c_const(qi(1)), Bq.den), x, lo, hi);
        Q al;
        if (c_const_value(c_div(integ_diff_poly(integ_diff_poly(q, x), x), c_const(qi(2))), &al) && q_sign(al) > 0)
        {                                             /* the hyperbola: q keeps its sign between the ends */
            Q a, b;
            if (!c_const_value(lo, &a) || !c_const_value(hi, &b)) nm_fail("ends of a hyperbola's area must be rational here");
            if (q_cmp(a, b) > 0) { Q t = a; a = b; b = t; }
            Q qa, qb;
            c_const_value(integ_subst_poly(q, x, c_const(a)), &qa); c_const_value(integ_subst_poly(q, x, c_const(b)), &qb);
            if (q_sign(qa) < 0 || q_sign(qb) < 0) nm_fail("the quantity under the root is negative at an end");
            Q mid = q_div(q_add(a, b), qi(2)), qm;
            c_const_value(integ_subst_poly(q, x, c_const(mid)), &qm);
            int inside = 0;
            if (q_sign(qa) == 0 || q_sign(qb) == 0) {  /* an end at a root: the other root, -be/al minus it */
                Q al2, be2, r;
                c_const_value(c_div(integ_diff_poly(integ_diff_poly(q, x), x), c_const(qi(2))), &al2);
                c_const_value(integ_subst_poly(integ_diff_poly(q, x), x, c_const(qi(0))), &be2);
                r = q_sub(q_neg(q_div(be2, al2)), q_sign(qa) == 0 ? a : b);
                inside = q_cmp(r, a) > 0 && q_cmp(r, b) < 0;
            } else inside = elim_has_root_closed(q, x, a, b);
            if (q_cmp(a, b) && (q_sign(qm) <= 0 || inside))
                nm_fail("the quantity under the root is negative inside the interval");
        }
        *out = vc(c_sub(sqrt_integral_value(I, q, hi, conic_constant), sqrt_integral_value(I, q, lo, conic_constant)));
        return 1;
    }
    if (!I.n) { *out = vr(I.rational); return 1; }
    Val r = vq(qi(0)); r.kind = V_AREA; r.area = I;
    *out = r;
    return 1;
}

static Val eval(Node *n) {
    switch (n->kind) {
    case N_MAT: {
        Mat m;
        const char *conditions = NULL;
        if (n->op == 'v') {                              /* a column */
            m = mat_new(n->nargs, 1);
            for (int i = 0; i < n->nargs; i++) {
                Val cell = eval(n->args[i]); m.a[i] = entry_of(cell);
                conditions = combine_verdicts(conditions,cell.verdict);
            }
        } else {
            int cols = n->args[0]->nargs;
            m = mat_new(n->nargs, cols);
            for (int i = 0; i < n->nargs; i++) {
                if (n->args[i]->nargs != cols) nm_fail("every row of a matrix must have %d entries", cols);
                for (int j = 0; j < cols; j++) {
                    Val cell = eval(n->args[i]->args[j]); m.a[i * cols + j] = entry_of(cell);
                    conditions = combine_verdicts(conditions,cell.verdict);
                }
            }
        }
        Val v = vmat(m); v.verdict = conditions; return v;
    }
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
    case N_BIN: {
        Val left = eval(n->a);
        Val right = eval(n->b);
        return arith(n->op, left, right);
    }
    case N_SQRT: return power_val(eval(n->a), vq(q_make(z_from_i64(1), z_from_i64(2))));
    case N_APPLY: {
        Binding *b = lookup(n->s, n->len), summed;
        if(b && b->v.kind==V_APART) {
            summed=*b;summed.v=vr(apart_sum(b->v.apart));summed.v.verdict=b->v.verdict;b=&summed;
        }
        if (b && b->v.kind == V_FUNC && !param_find(n->s, n->len)) return call_rule(b, n);
        if (!b && !param_find(n->s, n->len)) { Val bv; if (builtin(n->s, n->len, n, &bv)) return bv; }
        if (n->nargs > 1) nm_fail("%.*s is not a rule of several arguments", (int)n->len, n->s);
        if (b && b->v.kind == V_REC) {
            Val arg = eval(n->a);
            if (arg.kind == V_RAT || (arg.kind == V_POLY && c_has_plain(arg.c))) longjmp(need_series, 1);
            if (!recipe_is_equation(b)) {                    /* a definition by an expression: put the number in */
                if (n->primes) nm_fail("the fluxion of %s at a number: take d/d%s first, then put the number in", b->name, b->v.var);
                int saved_base = frame_base, saved_depth = pdepth;
                frame_base = pdepth;
                Param *p = param_push(b->v.var, strlen(b->v.var));
                p->hasv = 1; p->v = arg;
                Val r = eval(parse_text(b->v.src));
                pdepth = saved_depth; frame_base = saved_base;
                return r;
            }
            if (arg.kind == V_RAT || (arg.kind == V_POLY && c_has_plain(arg.c))) longjmp(need_series, 1);
            return series_value_d(b, arg, n->primes);
        }
        if (b && b->v.kind == V_AREA) {
            if(n->primes) nm_fail("take d/d%s first, then substitute into the conic area",letter_name(b->v.area.var));
            return vc(integ_value(b->v.area,as_c(eval(n->a)),conic_constant));
        }
        if (b && b->v.kind == V_RAT) {
            int li = rational_letter(b->v.rat);
            if (li >= 0) {
                Val arg = eval(n->a);
                if (!is_exact(arg)) nm_fail("a rational function currently takes an exact argument");
                if (n->primes) nm_fail("take d/d%s first, then substitute into the rational function",letter_name(li));
                R x = as_r(arg);
                R num = substitute_r(b->v.rat.num,li,x), den = substitute_r(b->v.rat.den,li,x);
                return vr(r_div(num,den));
            }
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
                    if (!q_is_int(ct_e(&term, li)) || !z_fits_i64(ct_e(&term, li).num, &e)) nm_fail("a fractional power at a number is not in this version");
                    acc = q_add(acc, q_mul(term.k, q_pow(arg.q, e)));
                }
                return vq(acc);
            }
        }
        Val left = name_val(n->s, n->len, n->primes);
        Val right = eval(n->a);
        return arith('*', left, right);
    }
    case N_DERIV: case N_INTEG: {
        if (n->kind == N_INTEG && !sqo.on) { Val rv; if (root_integral(n, &rv)) return rv; }
        Val v = eval(n->a);
        if(v.kind==V_APART) v=vr(apart_sum(v.apart));
        if(v.kind==V_AREA) {
            if(n->kind!=N_DERIV) nm_fail("integrating a conic area again comes later");
            if(v.area.root) nm_fail("the moment of an area under a root comes later; it is the integrand");
            if(n->len && letter_index(n->s,n->len)!=v.area.var) nm_fail("a conic area currently supports its own integration letter only");
            return vr(integ_derivative(v.area));
        }
        if(n->kind==N_INTEG && is_exact(v)) {
            int li=n->len?letter_index(n->s,n->len):rational_letter(as_r(v));
            int finite=li>=0;
            if(v.kind==V_POLY && li>=0) for(int t=0;t<v.c.nt;t++)
                for(int l=0;l<letter_count();l++) if(!q_is_int(ct_e(&v.c.t[t], l))) finite=0;
            if(finite) {
                R f=as_r(v);
                C lo=c_zero(),hi=c_zero();
                if(n->b) {
                    Val lv=eval(n->b),hv=eval(n->c);
                    if(!is_exact(lv) || !is_exact(hv)) nm_fail("definite integration needs exact real endpoints, not approximate values");
                    lo=as_c(lv);hi=as_c(hv);integral_no_poles(f,li,lo,hi);
                }
                Integral a=integ_rational(f,li);
                if(n->b) return vc(c_sub(integ_value(a,hi,conic_constant),integ_value(a,lo,conic_constant)));
                if(!a.n) return vr(a.rational);
                Val out=vq(qi(0));out.kind=V_AREA;out.area=a;out.verdict=v.verdict;return out;
            }
        }
        if(n->kind==N_INTEG && n->b) nm_fail("definite integration currently needs a rational function");
        if (v.kind == V_RAT && n->kind == N_DERIV) {
            int li = n->len ? letter_index(n->s,n->len) : rational_letter(v.rat);
            if (li < 0) longjmp(need_series,1);
            C dn = derivative_c(v.rat.num,li), dd = derivative_c(v.rat.den,li);
            return vr(r_make(c_sub(c_mul(dn,v.rat.den),c_mul(v.rat.num,dd)),c_pow_int(v.rat.den,2)));
        }
        if (v.kind == V_RAT && !c_is_monomial(v.rat.den)) longjmp(need_series, 1);
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
            Q e = ct_e(&term, li);
            if (n->kind == N_DERIV) {
                if (q_sign(e) == 0) continue;
                term.k = q_mul(term.k, e); ct_set(&term, li, q_sub(e, qi(1)));
            } else {                                       /* Newton's first rule: a x^m gives a x^(m+1)/(m+1) */
                if (q_cmp(e, qi(-1)) == 0) nm_fail("the fluent of 1/%s is a logarithm: write it as a root of an equation", letter_name(li));
                ct_set(&term, li, q_add(e, qi(1))); term.k = q_div(term.k, ct_e(&term, li));
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
    int has_sub; Dual sub;      /* the letter stands for this series (a definition used at another argument) */
    int has_shift; C shift;     /* the letter stands for shift + x (an equation moved to another point) */
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
static Ser resolve_ex(Node *root, const char *var, int n, const C *shift, const C *starts);
#define resolve(root, var, n) resolve_ex(root, var, n, NULL, NULL)
static int recipe_is_equation(Binding *b);
static Ser moved_series(Binding *b, C a, int n);
static CBall as_cball(Val v);

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
        if (!q_is_int(ct_e(&term, vi)) || !z_fits_i64(ct_e(&term, vi).num, &e) || e < 0)
            nm_fail("a negative or fractional power of %s: use the parallelogram (root of ... for y)", var);
        if (e >= n) continue;
        ct_set(&term, vi, qi(0));
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

/* a built-in at a number inside a series: if its arguments need the series letter, it has no series */
static int try_builtin(Node *n, Val *out) {
    jmp_buf saved; memcpy(saved, need_series, sizeof saved);
    volatile int ok = 0;
    if (!setjmp(need_series)) ok = builtin(n->s, n->len, n, out);
    memcpy(need_series, saved, sizeof saved);
    return ok;
}

static int recipe_is_equation(Binding *b) {
    Node *r = parse_text(b->v.src);
    return r->kind == N_ROOT && r->nconds;
}

/* the value of f^(j) at a, kept as a letter: exp(1), exp'(1) ... */
typedef struct { Binding *b; C a; int j; } MovedStart;

static Ball moved_start_value(void *data, int64_t prec) {
    MovedStart *m = data;
    int64_t saved = work_prec;
    work_prec = prec;
    Val av; Q k;
    if (c_const_value(m->a, &k)) av = vq(k); else { memset(&av, 0, sizeof av); av.kind = V_POLY; av.c = m->a; }
    CBall z = as_cball(series_value_d(m->b, av, m->j));
    work_prec = saved;
    if (z.im.m.s || z.im.r.s) nm_fail("a moved series at a complex point comes later");
    return z.re;
}

/* f(a + x) for f known by its equation: the same equation with x put for a + x, and its starts the values of
 * f, f', ... at a, each a letter known by its value */
static Ser moved_series(Binding *b, C a, int n) {
    Node *root = parse_text(b->v.src);
    const char *u = root->conds[0].s; size_t ul = root->conds[0].len;
    int r = max_primes(root->a, u, ul), r2 = max_primes(root->b, u, ul);
    if (r2 > r) r = r2;
    int need = r ? r : 1;
    if (c_has_plain(a)) nm_fail("moving %s to a point with letters comes later", b->name);
    C *starts = arena_alloc((size_t)need * sizeof(C));
    char *as = c_to_str(a);
    for (int j = 0; j < need; j++) {
        char *disp = arena_alloc(strlen(b->name) + strlen(as) + 16);
        { char *o = disp + sprintf(disp, "%s", b->name); for (int q = 0; q < j; q++) *o++ = 39; sprintf(o, "(%s)", as); }
        MovedStart *ms = perm_alloc(sizeof *ms);
        ms->b = b; ms->a = c_persist(a); ms->j = j;
        starts[j] = c_named(disp, moved_start_value, ms);
    }
    return resolve_ex(root, b->v.var, n, &a, starts);
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
        if (pp->v.kind == V_RAT) return dconst(s_div(s_from_c(pp->v.rat.num, cx->var, cx->n), s_from_c(pp->v.rat.den, cx->var, cx->n)));
        nm_fail("%.*s is an approximate number; it cannot enter a series", (int)len, s);
    }
    if (same(s, len, cx->var, strlen(cx->var))) {
        if (cx->has_sub) {
            Dual d = cx->sub;
            d.a = s_trunc(d.a, cx->n);
            if (d.hasb) d.b = s_trunc(d.b, cx->n);
            return d;
        }
        if (cx->has_shift) return dconst(s_add(s_var(cx->n), s_const(cx->shift, cx->n)));
        return dconst(s_var(cx->n));
    }
    Binding *b = lookup(s, len);
    if (!b && len == 1 && s[0] == 'i') return dconst(s_const(c_imag_unit(), cx->n));
    if (!b) return dconst(s_const(c_letter(letter_index(s, len)), cx->n));    /* a given letter: a, b, ... */
    switch (b->v.kind) {
    case V_Q: return dconst(s_const(c_const(b->v.q), cx->n));
    case V_POLY: return dconst(s_from_c(b->v.c, cx->var, cx->n));
    case V_RAT: return dconst(s_div(s_from_c(b->v.rat.num, cx->var, cx->n), s_from_c(b->v.rat.den, cx->var, cx->n)));
    case V_APART: {
        R r=apart_sum(b->v.apart);
        return dconst(s_div(s_from_c(r.num,cx->var,cx->n),s_from_c(r.den,cx->var,cx->n)));
    }
    case V_AREA: nm_fail("a stored conic area has no formal series in this version; expand integral(f,x) directly");
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
        if (v.kind == V_RAT) return dconst(s_div(s_from_c(v.rat.num, cx->var, cx->n), s_from_c(v.rat.den, cx->var, cx->n)));
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
        Binding *b = lookup(n->s, n->len), summed;
        if(b && b->v.kind==V_APART) {
            summed=*b;summed.v=vr(apart_sum(b->v.apart));summed.v.verdict=b->v.verdict;b=&summed;
        }
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
        if (!b && !param_find(n->s, n->len)) {
            Val bv;
            if (try_builtin(n, &bv)) {
                if(bv.kind==V_APART) bv=vr(apart_sum(bv.apart));
                if (bv.kind == V_Q) return dconst(s_const(c_const(bv.q), cx->n));
                if (bv.kind == V_POLY) return dconst(s_from_c(bv.c, cx->var, cx->n));
                if (bv.kind == V_RAT) return dconst(s_div(s_from_c(bv.rat.num, cx->var, cx->n), s_from_c(bv.rat.den, cx->var, cx->n)));
                nm_fail("%.*s gives no series", (int)n->len, n->s);
            }
        }
        if (n->nargs > 1) nm_fail("%.*s is not a rule of several arguments", (int)n->len, n->s);
        if (b && b->v.kind == V_RAT) {
            int li = rational_letter(b->v.rat);
            if (li >= 0) {
                if (n->primes) nm_fail("take d/d%s first, then substitute into the rational function",letter_name(li));
                Dual arg = sev(n->a,cx), values[2];
                C parts[2] = {b->v.rat.num,b->v.rat.den};
                for (int j = 0; j < 2; j++) {
                    UP p = up_from(parts[j],li);
                    values[j] = dconst(s_const(c_zero(),cx->n));
                    for (int i = p.deg; i >= 0; i--)
                        values[j] = dadd(dmul(values[j],arg),dconst(s_from_c(p.c[i],cx->var,cx->n)),1);
                }
                return ddiv(values[0],values[1]);
            }
        }
        if (!b || b->v.kind != V_REC) {
            Dual left = sname(n->s, n->len, n->primes, cx);
            Dual right = sev(n->a, cx);
            return dmul(left, right);
        }
        Dual g = sev(n->a, cx);
        if (!recipe_is_equation(b)) {                        /* a definition by an expression: put the series in */
            SCtx sc; memset(&sc, 0, sizeof sc);
            sc.var = b->v.var; sc.n = cx->n; sc.has_sub = 1; sc.sub = g;
            return sev(parse_text(b->v.src), &sc);
        }
        if (g.a.n && !c_is_zero(g.a.c[0])) {                /* f(a + ...): the equation moved to a */
            C a0 = g.a.c[0];
            Ser f = moved_series(b, a0, g.a.n);
            Dual h = g;
            h.a = s_sub(g.a, s_const(a0, g.a.n));
            Dual d;
            d.a = s_compose(f, h.a);
            d.hasb = g.hasb;
            if (g.hasb) d.b = s_mul(s_compose(s_deriv(f), s_trunc(h.a, f.n - 1)), g.b);
            return d;
        }
        Ser f = recipe_series(b, b->v.var, g.a.n);                /* f(g), and its moment f'(g) g' */
        Dual d;
        d.a = s_compose(f, g.a);
        d.hasb = g.hasb;
        if (g.hasb) d.b = s_mul(s_compose(s_deriv(f), s_trunc(g.a, f.n - 1)), g.b);
        return d;
    }
    case N_DERIV: case N_INTEG: {
        if(n->kind==N_INTEG && n->b) nm_fail("definite integrals take places, not a formal series order");
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
static Ser resolve_ex(Node *root, const char *var, int n, const C *shift, const C *starts) {
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
    if (starts) {                                       /* the equation moved: its starts are given values */
        for (int j = 0; j < need; j++) {
            Q fact = qi(1);
            for (int f = 2; f <= j; f++) fact = q_mul(fact, qi(f));
            c[j] = c_scale(starts[j], q_div(qi(1), fact));
            given[j] = 1;
        }
    }
    for (int i = 0; i < root->nconds && !starts; i++) {
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
    if (shift) { cx.has_shift = 1; cx.shift = *shift; }
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
typedef struct { const char *s[64]; size_t len[64]; int n; } Letters;

static void add_letter(Letters *L, const char *s, size_t len) {
    for (int i = 0; i < L->n; i++) if (same(L->s[i], L->len[i], s, len)) return;
    if (L->n < 64) { L->s[L->n] = s; L->len[L->n] = len; L->n++; }
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
        else if (b && b->v.kind == V_RAT) {
            for (int l = 0; l < letter_count(); l++) if (c_uses(b->v.rat.num, l) || c_uses(b->v.rat.den, l))
                add_letter(L, letter_name(l), strlen(letter_name(l)));
        } else if (b && b->v.kind == V_POLY) {
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
typedef TermRule Rule;

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

/* the value of f, or of its j-th fluxion, at a number: real or complex, exact or approximate */
static Val series_value_d(Binding *b, Val arg, int deriv) {
    CBall A = as_cball(arg);
    int real = !A.im.m.s && !A.im.r.s;
    Q a = arg.kind == V_Q ? arg.q : q_from_z(z_zero());
    Node *root = parse_text(b->v.src);
    if (root->kind != N_ROOT || !root->nconds)
        nm_fail("the value of %s needs its definition as an equation with a start (root of ..., y(0) = ...)", b->name);
    const char *u = root->conds[0].s; size_t ul = root->conds[0].len;
    int r = max_primes(root->a, u, ul), r2 = max_primes(root->b, u, ul);
    if (r2 > r) r = r2;
    Rule R = term_rule(root, b->v.var, r);
    if (R.D.deg < 0) nm_fail("the equation fixes no term of %s", b->name);
    Q rr = arg.kind == V_Q ? qabs(a) : cb_abs_upper(A);
    /* how fast the terms shrink in the end: sum over t of lim |N_t / D| r^t must be below 1 */
    Q ginf = qi(0);
    for (int t = 1; t <= R.T; t++) {
        if (R.N[t].deg > R.D.deg) nm_fail("the terms of %s grow too fast: the series has no value away from 0", b->name);
        if (R.N[t].deg == R.D.deg) ginf = q_add(ginf, q_mul(qabs(q_div(R.N[t].c[R.D.deg], R.D.c[R.D.deg])), q_pow(rr, t)));
    }
    if (q_cmp_one(ginf) >= 0)
        nm_fail("at %s the terms of %s shrink too slowly to bound the rest; use smaller arguments and exact relations, as Newton did with 0.1 and 0.2",
                arg.kind == V_Q ? q_to_str(a) : "this point", b->name);
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
    if (arg.kind == V_Q && q_sign(a) == 0) {                /* at 0: the term itself, times j! */
        while (ncoef <= deriv) {
            int n = ncoef;
            Q acc = qi(0);
            for (int t = 1; t <= R.T; t++) acc = q_add(acc, q_mul(p_at(R.N[t], qi(n)), c[n - t]));
            int m = n - R.s;
            if (m >= 0 && m <= R.h.deg) acc = q_sub(acc, R.h.c[m]);
            if (ncoef == cap) nm_fail("too many terms");
            c[n] = q_div(acc, p_at(R.D, qi(n)));
            ncoef++;
        }
        Q f = c[deriv];
        for (int u = 2; u <= deriv; u++) f = q_mul(f, qi(u));
        return vq(f);
    }
    CBall v = rule_value(&R, c, ncoef, start, A, deriv, work_prec - GUARD_DIGITS, work_prec);
    if (real) return vball(v.re);
    return vcb(v);
}

/* The same equations as lib/prelude.nm, through the existing resolver and
 * certified term rules. Private bindings keep mathematical constants stable
 * when a user gives atan, log1p, x or y another meaning. */
static Ball conic_series(int kind,Ball x,int64_t prec) {
    Binding b={0}; b.name=kind==AREA_LOG?"log1p":"atan";
    b.v.kind=V_REC;b.v.var="x";
    b.v.src=kind==AREA_LOG?"root of (1 + x) y' = 1, y(0) = 0":"root of (1 + x^2) y' = 1, y(0) = 0";
    Binding *save_names=names;int save_base=frame_base;int64_t save_prec=work_prec;
    jmp_buf saved_error;memcpy(saved_error,nm_on_error,sizeof saved_error);
    if(setjmp(nm_on_error)) {
        char err[sizeof nm_error_msg];memcpy(err,nm_error_msg,sizeof err);
        names=save_names;frame_base=save_base;work_prec=save_prec;
        memcpy(nm_on_error,saved_error,sizeof saved_error);nm_fail("%s",err);
    }
    names=NULL;frame_base=pdepth;work_prec=prec+GUARD_DIGITS;
    Ball r=as_ball(series_value_d(&b,vball(x),0));
    names=save_names;frame_base=save_base;work_prec=save_prec;
    memcpy(nm_on_error,saved_error,sizeof saved_error);
    return r;
}
static Q conic_end(Ball b,int upper) {
    Z m=upper?z_add(b.m,b.r):z_sub(b.m,b.r);
    return b.e>=0?q_from_z(z_mul_pow10(m,b.e)):q_make(m,z_pow10(-b.e));
}
static int number_sign(C x) {
    if(c_is_zero(x)) return 0;
    if(c_has_plain(x)) nm_fail("a conic value requires real numbers, without parameters");
    Ball b=ball_of_c(x,work_prec+10);
    if(q_sign(conic_end(b,0))>0) return 1;
    if(q_sign(conic_end(b,1))<0) return -1;
    nm_fail("the sign of a conic argument is not separated at this precision");
}
static Ball conic_pi(void *data,int64_t prec) {
    /* atan(1)=2 atan(1/3)+atan(1/7); all series arguments are inside the disc. */
    int64_t p=prec+12;
    Ball a=conic_series(AREA_ATAN,b_from_q(q_div(qi(1),qi(3)),p),p);
    Ball b=conic_series(AREA_ATAN,b_from_q(q_div(qi(1),qi(7)),p),p);
    return b_add(b_mul(b_from_q(qi(8),p),a,p),b_mul(b_from_q(qi(4),p),b,p),prec);
}
typedef struct { int kind; C arg; } ConicNumber;
static Ball conic_number(void *data,int64_t prec) {
    ConicNumber *c=data;int64_t p=prec+12;
    Ball x=ball_of_c(c->arg,p),one=b_from_q(qi(1),p),two=b_from_q(qi(2),p);
    if(c->kind==AREA_LOG) {
        int shift=0;
        while(q_cmp(conic_end(x,1),q_div(qi(3),qi(2)))>0) {
            x=b_div(x,two,p);if(++shift>100000) nm_fail("logarithm range reduction limit");
        }
        while(q_cmp(conic_end(x,0),q_div(qi(3),qi(4)))<0) {
            x=b_mul(x,two,p);if(--shift < -100000) nm_fail("logarithm range reduction limit");
        }
        Ball r=conic_series(AREA_LOG,b_sub(x,one,p),p);
        if(shift) {
            /* log(2)=log1p(1/2)-log1p(-1/4). */
            Ball a=conic_series(AREA_LOG,b_from_q(q_div(qi(1),qi(2)),p),p);
            Ball b=conic_series(AREA_LOG,b_from_q(q_div(qi(-1),qi(4)),p),p);
            r=b_add(r,b_mul(b_from_q(qi(shift),p),b_sub(a,b,p),p),p);
        }
        return r;
    }
    /* Positive arguments only. Invert above 2; subtract atan(1/2) above 1/2.
     * The resulting absolute argument is at most 3/4. */
    int invert=q_cmp(conic_end(x,0),qi(2))>0;
    if(invert) x=b_div(one,x,p);
    Ball r;
    if(q_cmp(conic_end(x,1),q_div(qi(1),qi(2)))>0) {
        Ball half=b_from_q(q_div(qi(1),qi(2)),p);
        Ball y=b_div(b_sub(x,half,p),b_add(one,b_mul(x,half,p),p),p);
        r=b_add(conic_series(AREA_ATAN,half,p),conic_series(AREA_ATAN,y,p),p);
    } else r=conic_series(AREA_ATAN,x,p);
    if(invert) r=b_sub(b_div(conic_pi(NULL,p),two,p),r,p);
    return r;
}
static C conic_constant(int kind,C arg) {
    int sign=number_sign(arg);
    if(!sign) {
        if(kind==AREA_LOG) nm_fail("logarithm at zero: pole");
        return c_zero();
    }
    if(sign<0) arg=c_neg(arg);
    if(kind==AREA_LOG && c_equal(arg,c_const(qi(1)))) return c_zero();
    if(kind==AREA_ATAN) {
        Q sq;C square=c_mul(arg,arg);
        if(c_const_value(square,&sq)) {
            Q k=qi(0);
            if(q_cmp(sq,qi(1))==0) k=q_div(qi(sign),qi(4));
            if(q_cmp(sq,qi(3))==0) k=q_div(qi(sign),qi(3));
            if(q_cmp(sq,q_div(qi(1),qi(3)))==0) k=q_div(qi(sign),qi(6));
            if(q_sign(k)) return c_scale(c_named("pi",conic_pi,NULL),k);
        }
    }
    ConicNumber *n=perm_alloc(sizeof *n);n->kind=kind;n->arg=c_persist(arg);
    char *a=c_to_str(arg),*label=arena_alloc(strlen(a)+10);
    sprintf(label,"%s(%s)",kind==AREA_LOG?"log":"atan",a);
    C r=c_named(label,conic_number,n);
    return kind==AREA_ATAN && sign<0?c_neg(r):r;
}
static void integral_no_poles(R f,int v,C lo,C hi) {
    /* Exact endpoint substitution detects poles at algebraic endpoints. A
     * Sturm count over enclosing rational endpoints then certifies the interval. */
    number_sign(lo);number_sign(hi);
    if(c_is_zero(integ_subst_poly(f.den,v,lo)) || c_is_zero(integ_subst_poly(f.den,v,hi)))
        nm_fail("pole in the closed integration interval");
    if(c_equal(lo,hi)) return;
    Q exact_lo,exact_hi;
    if(c_const_value(lo,&exact_lo) && c_const_value(hi,&exact_hi)) {
        if(elim_has_root_closed(f.den,v,exact_lo,exact_hi)) nm_fail("pole in the closed integration interval");
        return;
    }
    int order=number_sign(c_sub(hi,lo));
    if(order<0) { C t=lo;lo=hi;hi=t; }
    Ball lb=ball_of_c(lo,work_prec+10),hb=ball_of_c(hi,work_prec+10);
    Q l=conic_end(lb,0),h=conic_end(hb,1);
    if(!elim_has_root_closed(f.den,v,l,h)) return;
    Q inner_l=conic_end(lb,1),inner_h=conic_end(hb,0);
    if(q_cmp(inner_l,inner_h)<=0 && elim_has_root_closed(f.den,v,inner_l,inner_h))
        nm_fail("pole in the closed integration interval");
    nm_fail("a pole is too close to an endpoint to separate at this precision");
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
    if (c->op == '=' && is_exact(a) && is_exact(b)) return r_equal(as_r(a), as_r(b));
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
    case V_AREA: {
        body=integ_to_str(v.area); snprintf(out,4096,"[exact]");
        C *vals=arena_alloc((size_t)(2*v.area.n+2)*sizeof(C)); int nv=0;          /* roots r_k met in the area */
        for(int i=0;i<v.area.n;i++) { vals[nv++]=v.area.term[i].poly; vals[nv++]=c_div(v.area.term[i].coef.num,v.area.term[i].coef.den); }
        char *leg=arena_alloc(4096),*lo=leg; *leg=0;
        root_legend(vals,nv,&lo);
        if(*leg) { char *o2=arena_alloc(strlen(out)+strlen(leg)+2); sprintf(o2,"%s%s",out,leg); snprintf(out,4096,"%s",o2); }
        break;
    }
    case V_APART: body=apart_to_str(v.apart); snprintf(out,4096,"[exact]"); break;
    case V_RAT: body = r_to_str(v.rat); snprintf(out, 4096, "[exact]"); break;
    case V_POLY: {
        C re, im;
        int lone = v.c.nt == 1 && q_is_int(v.c.t[0].k) && z_is_one(v.c.t[0].k.num);
        int li = -1, nl = 0;
        for (int l = 0; l < letter_count(); l++) if (c_uses(v.c, l)) { li = l; nl++; }
        if (asked && lone && nl == 1 && letter_is_surd(li) && surd_root(li) && q_cmp_one(ct_e(&v.c.t[0], li)) == 0) {
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
        nm_fact("equation", "%s", nm_json_str(root_equation_str(r)));
        if (r->certified && r->w)                    /* the sign change was seen across this interval */
            nm_fact("interval", "[%s,%s]", nm_json_str(fixed_str(z_sub(r->X, z_from_i64(r->w)), r->D)),
                    nm_json_str(fixed_str(z_add(r->X, z_from_i64(r->w)), r->D)));
        else if (r->certified) nm_fact("exact_root", "true");
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
    case V_MAT: body = mat_to_str(v.m); if (!v.verdict) v.verdict = mat_verdict(v.m, mat_new(0,0), -1); snprintf(out, 4096, "[exact]"); break;
    case V_TEXT: return (char *)v.text;
    case V_CBALL: {
        int64_t kr = b_guaranteed_places(v.z.re, places), ki = b_guaranteed_places(v.z.im, places);
        int64_t k = kr < ki ? kr : ki;
        if (k < 0) nm_fail("not even the units place is guaranteed");
        char *sr = fixed_str(b_round_to_places(v.z.re, k), k), *si = fixed_str(b_round_to_places(v.z.im, k), k);
        body = arena_alloc(strlen(sr) + strlen(si) + 8);
        if (si[0] == '-') sprintf(body, "%s - %si", sr, si + 1); else sprintf(body, "%s + %si", sr, si);
        if (k == places) snprintf(out, 4096, "[bounded: %lld places guaranteed]", (long long)k);
        else snprintf(out, 4096, "[bounded: only %lld of %lld places guaranteed]", (long long)k, (long long)places);
        break;
    }
    default: body = "?";
    }
    const char *verdict = v.verdict ? v.verdict : out;
    if (v.verdict && strcmp(out,"[exact]")) {
        char *qualified = arena_alloc(strlen(out) + strlen(v.verdict) + 24);
        sprintf(qualified,"%s [subject to %s]",out,v.verdict);
        verdict = qualified;
    }
    char *line = arena_alloc(strlen(body) + strlen(verdict) + 4);
    sprintf(line, "%s  %s", body, verdict);
    return line;
}

static Val persist_val(Val v) {
    if (v.verdict) v.verdict = dup_perm(v.verdict, strlen(v.verdict));
    switch (v.kind) {
    case V_Q: v.q = q_persist(v.q); break;
    case V_POLY: v.c = c_persist(v.c); break;
    case V_AREA: v.area=integ_persist(v.area); break;
    case V_APART: v.apart=apart_persist(v.apart); break;
    case V_RAT: v.rat = r_persist(v.rat); break;
    case V_BALL: v.ball.m = z_persist(v.ball.m); v.ball.r = z_persist(v.ball.r); break;
    case V_CBALL:
        v.z.re.m = z_persist(v.z.re.m); v.z.re.r = z_persist(v.z.re.r);
        v.z.im.m = z_persist(v.z.im.m); v.z.im.r = z_persist(v.z.im.r);
        break;
    case V_MAT: {
        Mat m; m.r = v.m.r; m.c = v.m.c;
        m.a = perm_alloc((size_t)(m.r * m.c > 0 ? m.r * m.c : 1) * sizeof(R));
        for (int i = 0; i < m.r * m.c; i++) m.a[i] = r_persist(v.m.a[i]);
        v.m = m;
        break;
    }
    case V_TEXT: { size_t L = strlen(v.text); char *t = perm_alloc(L + 1); memcpy(t, v.text, L + 1); v.text = t; break; }
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

char *nm_run(const char *line, int *failed);

extern const char *nm_lib_names[], *nm_lib_texts[];

/* use NAME: the library from NEWTONMATH_LIB if set, else the one compiled in, else NAME.nm in this directory */
static char *use_library(const char *name, size_t len) {
    const char *dir = getenv("NEWTONMATH_LIB");
    char path[1024];
    FILE *f = NULL;
    if (dir) { snprintf(path, sizeof path, "%s/%.*s.nm", dir, (int)len, name); f = fopen(path, "r"); }
    for (int i = 0; !f && nm_lib_names[i]; i++)
        if (strlen(nm_lib_names[i]) == len && !strncmp(nm_lib_names[i], name, len)) {
            snprintf(path, sizeof path, "%.*s (built in)", (int)len, name);
            f = fmemopen((void *)nm_lib_texts[i], strlen(nm_lib_texts[i]), "r");
        }
    if (!f) { snprintf(path, sizeof path, "%.*s.nm", (int)len, name); f = fopen(path, "r"); }
    if (!f) nm_fail("no library %.*s (built in or %s)", (int)len, name, path);
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
/* ---------------- systems: solve ... for x, y  and  eliminate y from ... ---------------- */

static C equation_c(void) {                            /* A = B  as  A - B, a sum of letters */
    Node *a = expr(), *b = NULL;
    if (at_op('=')) { pos++; b = expr(); }
    jmp_buf saved; memcpy(saved, need_series, sizeof saved);
    if (setjmp(need_series)) { memcpy(need_series, saved, sizeof saved); nm_fail("the equations of a system are sums of letters and numbers"); }
    Val va = eval(a), vb = b ? eval(b) : vq(qi(0));
    memcpy(need_series, saved, sizeof saved);
    Val d = arith('-', va, vb);
    if (d.kind == V_Q) return c_const(d.q);
    if (d.kind == V_POLY) return d.c;
    nm_fail("the equations of a system must be exact");
}

static int equations(C *eqs, int max) {
    int ne = 0;
    for (;;) {
        if (ne == max) nm_fail("at most %d equations", max);
        eqs[ne++] = equation_c();
        if (!at_op(',')) break;
        pos++;
    }
    return ne;
}

/* the roots r_1, r_2 ... met in the answer, each with its equation and certified places */
static void root_legend(C *vals, int nv, char **o) {
    for (int l = 0; l < letter_count(); l++) {
        if (!letter_is_surd(l) || !surd_root(l) || strchr(letter_name(l), '(')) continue;   /* sqrt(2) needs no legend */
        int used = 0;
        for (int j = 0; j < nv && !used; j++) used = c_uses(vals[j], l);
        if (!used) continue;
        Val rv; memset(&rv, 0, sizeof rv); rv.kind = V_ROOT; rv.root = surd_root(l);
        *o += sprintf(*o, "\n  where %s = %s", letter_name(l), show(rv, 20, 1));
    }
}

static char *solve_linear_letters(C *eqs, int ne, int *vars, int nv) {
    Mat a = mat_new(ne,nv), rhs = mat_new(ne,1);
    int parameters = 0;
    for (int i = 0; i < ne; i++) for (int t = 0; t < eqs[i].nt; t++) {
        CT term = eqs[i].t[t];
        int unknown = -1;
        for (int j = 0; j < nv; j++) if (q_sign(ct_e(&term, vars[j]))) {
            if (unknown >= 0 || q_cmp_one(ct_e(&term, vars[j])) != 0) return NULL;
            unknown = j; ct_set(&term, vars[j], qi(0));
        }
        C coefficient = {1,&term};
        if (c_has_plain(coefficient)) parameters = 1;
        R r = r_from_c(coefficient);
        if (unknown < 0) rhs.a[i] = r_sub(rhs.a[i],r);
        else a.a[i * nv + unknown] = r_add(a.a[i * nv + unknown],r);
    }
    if (!parameters) return NULL;
    Mat ns, x = mat_solve(a,rhs,&ns);
    Mat formulas = mat_new(nv,1 + ns.c);
    for (int i = 0; i < nv; i++) {
        formulas.a[i * formulas.c] = x.a[i];
        for (int j = 0; j < ns.c; j++) formulas.a[i * formulas.c + j + 1] = ns.a[i * ns.c + j];
    }
    char *verdict = mat_verdict(a,formulas,nv - ns.c), *basis = mat_to_str(ns);
    size_t cap = strlen(verdict) + strlen(basis) + 128;
    for (int i = 0; i < nv; i++) cap += strlen(letter_name(vars[i])) + strlen(r_to_str(x.a[i])) + 8;
    char *out = arena_alloc(cap), *o = out;
    for (int i = 0; i < nv; i++) o += sprintf(o,"%s%s = %s",i ? ", " : "",letter_name(vars[i]),r_to_str(x.a[i]));
    if (ns.c) o += sprintf(o,"; add any combination of the columns of %s",basis);
    sprintf(o,"  %s",verdict);
    return out;
}

static char *solve_stmt(void) {
    C eqs[16];
    int ne = equations(eqs, 16);
    expect_key("for");
    int vars[NM_MAXL], nv = 0;
    for (;;) {
        if (peek()->kind != T_NAME) nm_fail("expected an unknown after 'for'");
        if (nv == NM_MAXL) nm_fail("too many unknowns");
        vars[nv++] = letter_index(peek()->s, peek()->len);
        pos++;
        if (!at_op(',')) break;
        pos++;
    }
    if (peek()->kind != T_END) nm_fail("unexpected '%.*s'", (int)peek()->len, peek()->s);
    work_prec = DEFAULT_PLACES + GUARD_DIGITS;
    char *linear = solve_linear_letters(eqs,ne,vars,nv);
    if (linear) return linear;
    Solutions S = elim_solve(eqs, ne, vars, nv);
    size_t cap = 512;
    for (int k = 0; k < S.n; k++) for (int j = 0; j < nv; j++) cap += strlen(c_to_str(S.val[k][j])) + 400;
    char *txt = arena_alloc(cap), *o = txt;
    for (int k = 0; k < S.n; k++) {
        for (int j = 0; j < nv; j++) o += sprintf(o, "%s%s = %s", j ? ", " : "", letter_name(vars[j]), c_to_str(S.val[k][j]));
        root_legend(S.val[k], nv, &o);
        o += sprintf(o, "\n");
    }
    if (S.n == 0 && !S.curve && !S.unresolved) { sprintf(o, "no solution  [exact: the extermination ends in a number not 0]"); return txt; }
    if (S.n) o += sprintf(o, "[%d solution%s, each put back into the equations: exact", S.n, S.n == 1 ? "" : "s");
    else o += sprintf(o, "[no isolated solution");
    if (S.rejected) o += sprintf(o, "; %d brought in by the extermination and refused", S.rejected);
    if (S.unresolved) o += sprintf(o, "; and %d complex roots of factors above degree two, not shown", S.unresolved);
    if (S.curve) o += sprintf(o, "; some solutions form a curve, not points: not listed");
    sprintf(o, "]");
    return txt;
}

static C eliminate_expr(void) {
    if (peek()->kind != T_NAME) nm_fail("expected the letter to exterminate after 'eliminate'");
    int v = letter_index(peek()->s, peek()->len);
    pos++;
    expect_key("from");
    C eqs[16];
    int ne = equations(eqs, 16);
    if (ne != 2) nm_fail("eliminate takes two equations (Newton: from two equations one unknown is exterminated)");
    return elim_resultant(eqs[0], eqs[1], v);
}

static char *run_statement(const char *line, int *failed);

/* show STATEMENT: the same statement, with each step of the work written as it is done */
char *nm_run(const char *line, int *failed) {
    nm_show_work = 0;
    nm_account_reset();
    char *out = run_statement(line, failed);
    if (nm_show_work && !nm_work_lines && !*failed) nm_work("(no steps are written for this kind of statement yet)");
    if (nm_show_work) nm_account_print_facts();
    nm_show_work = 0;
    return out;
}

static char *run_statement(const char *line, int *failed) {
    *failed = 0;
    if (setjmp(nm_on_error)) {
        *failed = 1;
        char *m = arena_alloc(strlen(nm_error_msg) + 8);
        sprintf(m, "error: %s", nm_error_msg);
        return m;
    }
    integ_sign = number_sign;
    lex(line);
    pdepth = frame_base = calldepth = 0; nsum_skip = 0;
    stmt_id++;
    if (at_key("show")) { pos++; nm_show_work = 1; }
    if (peek()->kind == T_END) return NULL;
    if (at_key("use")) {
        pos++;
        if (peek()->kind != T_NAME) nm_fail("expected a library name after 'use'");
        Tok t = *peek();
        char *name = dup_perm(t.s, t.len);             /* the library's own lines replace this one */
        return use_library(name, t.len);
    }
    if (at_key("solve")) { pos++; return solve_stmt(); }
    const char *volatile let_name = NULL; volatile size_t let_len = 0;
    if (at_key("let")) {
        pos++;
        if (peek()->kind != T_NAME) nm_fail("expected a name after 'let'");
        let_name = peek()->s; let_len = peek()->len; pos++;
        if (at_op('(') || at_op('[')) return define_rule(let_name, let_len);
        expect_op('=');
    }
    if (at_key("eliminate")) {
        pos++;
        C r = eliminate_expr();
        if (peek()->kind != T_END) nm_fail("unexpected '%.*s'", (int)peek()->len, peek()->s);
        Val v; Q k;
        if (c_const_value(r, &k)) v = vq(k); else { memset(&v, 0, sizeof v); v.kind = V_POLY; v.c = r; }
        char *text = show(v, DEFAULT_PLACES, 0);
        if (!let_name) return text;
        bind(let_name, let_len, v);
        char *o = arena_alloc(strlen(text) + let_len + 4);
        sprintf(o, "%.*s = %s", (int)let_len, let_name, text);
        return o;
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
        if (to_var) longjmp(need_series, 1);
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
