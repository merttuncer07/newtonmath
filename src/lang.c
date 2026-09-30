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

static const char *KEYWORDS[] = {"let", "to", "places", "root", "of", "near", "sqrt", "integral", "use", NULL};

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
        if (strchr("+-*/^()=,", c)) { push(T_OP, src + i, 1, (char)c); i++; continue; }
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

enum { N_NUM, N_NAME, N_NEG, N_BIN, N_SQRT, N_ROOT, N_APPLY, N_DERIV, N_INTEG };

typedef struct Cond { const char *s; size_t len; int primes; struct Node *at, *val; } Cond;

typedef struct Node {
    int kind; char op;
    const char *s; size_t len; int primes;      /* a name, or the letter of d/dx and integral */
    struct Node *a, *b, *c;
    Cond *conds; int nconds;
} Node;

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
        if (at_op('(')) {                                              /* f(x), or y (x + 1) */
            pos++;
            Node *n = mk(N_APPLY);
            n->s = t->s; n->len = t->len; n->primes = t->primes;
            n->a = expr(); expect_op(')');
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

/* parse a stored definition without disturbing the statement being read */
static Node *parse_text(const char *src) {
    Tok *st = toks; int sn = ntok, sp = pos;
    lex(src);
    Node *e = expr();
    if (peek()->kind != T_END) nm_fail("internal: stored definition did not parse");
    toks = st; ntok = sn; pos = sp;
    return e;
}

/* ---------------- values ---------------- */

enum { V_Q, V_POLY, V_ROOT, V_BALL, V_REC, V_NUMREC };   /* recipes: a series, an approximate number */
typedef struct { int kind; Q q; Poly p; Root *root; Ball ball; const char *src, *var; } Val;

typedef struct Binding { char *name; Val v; struct Binding *next; } Binding;
static Binding *names;

static Binding *lookup(const char *s, size_t len) {
    for (Binding *b = names; b; b = b->next)
        if (strlen(b->name) == len && !strncmp(b->name, s, len)) return b;
    return NULL;
}

static int same(const char *a, size_t la, const char *b, size_t lb) { return la == lb && !strncmp(a, b, la); }

static char *dup_arena(const char *s, size_t len) { char *r = arena_alloc(len + 1); memcpy(r, s, len); r[len] = 0; return r; }
static char *dup_perm(const char *s, size_t len) { char *r = perm_alloc(len + 1); memcpy(r, s, len); r[len] = 0; return r; }

static Val vq(Q q) { Val v; memset(&v, 0, sizeof v); v.kind = V_Q; v.q = q; return v; }
static Val vpoly(Poly p) {
    if (p.deg <= 0) return vq(p.deg < 0 ? q_from_z(z_zero()) : p.c[0]);
    Val v; memset(&v, 0, sizeof v); v.kind = V_POLY; v.p = p; return v;
}
static Val vball(Ball b) { Val v; memset(&v, 0, sizeof v); v.kind = V_BALL; v.ball = b; return v; }

static int64_t work_prec;                               /* significant digits for balls in this statement */
static jmp_buf need_series;                             /* the exact pass meets something only a series can hold */

static Poly as_poly(Val v) { return v.kind == V_POLY ? v.p : p_const(v.q); }

static Ball as_ball(Val v) {
    switch (v.kind) {
    case V_Q: return b_from_q(v.q, work_prec);
    case V_ROOT: return b_from_root(v.root, work_prec);
    case V_BALL: return v.ball;
    default: nm_fail("an unknown letter cannot be mixed with approximate numbers");
    }
    return v.ball;
}

static int is_exact(Val v) { return v.kind == V_Q || v.kind == V_POLY; }
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

/* c^(1/n) for a rational c: exact when c is a perfect n-th power, otherwise a root of B*y^n - A */
static Val nth_root(Q c, int64_t n) {
    if (n < 2 || n > 1000) nm_fail("root index out of range");
    int neg = q_sign(c) < 0;
    if (neg && n % 2 == 0) nm_fail("an even root of a negative number is not real (complex numbers come later)");
    Q exact;
    if (q_root_exact(c, n, &exact)) return vq(exact);
    Z A = z_abs(c.num), B = c.den;
    int64_t D0 = 6;                                     /* start from the integer root, as by hand */
    Z scaled = z_iroot(z_div_round(z_mul_pow10(A, D0 * n), B), (unsigned)n);
    Q guess = q_make(neg ? z_neg(scaled) : scaled, z_pow10(D0));
    Poly p = p_pow(p_var("y"), (unsigned)n);
    p = p_sub(p_scale(p, q_from_z(B)), p_const(q_from_z(neg ? z_neg(A) : A)));
    Val v; memset(&v, 0, sizeof v);
    v.kind = V_ROOT;
    v.root = root_new(p, guess, D0);
    return v;
}

/* ---------------- the exact pass: numbers, rationals, polynomials ---------------- */

static Val eval(Node *n);
static Val series_value(Binding *b, Q a);

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
        if (den != 1 || num < 0) longjmp(need_series, 1);
        return vpoly(p_pow(base.p, (unsigned)num));
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
        Poly x = as_poly(a), y = as_poly(b);
        switch (op) {
        case '+': return vpoly(p_add(x, y));
        case '-': return vpoly(p_sub(x, y));
        case '*': return vpoly(p_mul(x, y));
        case '/':
            if (b.kind != V_Q) longjmp(need_series, 1);
            return vpoly(p_scale(x, q_div(qi(1), b.q)));
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
    Binding *b = lookup(s, len);
    if (b) {
        if (b->v.kind == V_REC) longjmp(need_series, 1);
        if (b->v.kind == V_NUMREC) return eval(parse_text(b->v.src));   /* carried again to the places now asked */
        return b->v;
    }
    return vpoly(p_var(dup_arena(s, len)));            /* an unknown letter */
}

static Val eval(Node *n) {
    switch (n->kind) {
    case N_NUM: return vq(parse_number(n->s, n->len));
    case N_NAME: return name_val(n->s, n->len, n->primes);
    case N_NEG: return arith('-', vq(qi(0)), eval(n->a));
    case N_BIN: return arith(n->op, eval(n->a), eval(n->b));
    case N_SQRT: return power_val(eval(n->a), vq(q_make(z_from_i64(1), z_from_i64(2))));
    case N_APPLY: {
        Binding *b = lookup(n->s, n->len);
        if (b && b->v.kind == V_REC) {
            Val arg = eval(n->a);
            if (arg.kind == V_POLY) longjmp(need_series, 1);
            if (arg.kind != V_Q) nm_fail("the value of %s needs an exact number as its argument in this version", b->name);
            return series_value(b, arg.q);
        }
        if (b && b->v.kind == V_POLY) {                     /* a polynomial at a number */
            Val arg = eval(n->a);
            if (arg.kind == V_Q) {
                Q acc = q_from_z(z_zero());
                for (int i = b->v.p.deg; i >= 0; i--) acc = q_add(q_mul(acc, arg.q), b->v.p.c[i]);
                return vq(acc);
            }
        }
        return arith('*', name_val(n->s, n->len, n->primes), eval(n->a));
    }
    case N_DERIV: case N_INTEG: {
        Val v = eval(n->a);
        if (!is_exact(v)) { if (n->kind == N_DERIV) return vq(qi(0)); longjmp(need_series, 1); }
        Poly p = as_poly(v);
        if (n->kind == N_INTEG && !n->len && !p.var) longjmp(need_series, 1);
        const char *var = n->len ? dup_arena(n->s, n->len) : p.var;
        if (p.var && strcmp(p.var, var)) {
            if (n->kind == N_DERIV) return vq(qi(0));
            longjmp(need_series, 1);
        }
        Poly r; r.var = var;
        if (n->kind == N_DERIV) {
            r.deg = p.deg - 1;
            r.c = arena_alloc((size_t)(p.deg > 0 ? p.deg : 1) * sizeof(Q));
            for (int i = 1; i <= p.deg; i++) r.c[i - 1] = q_mul(p.c[i], qi(i));
        } else {                                       /* Newton's first rule: a x^m gives a x^(m+1)/(m+1) */
            r.deg = p.deg + 1;
            r.c = arena_alloc((size_t)(p.deg + 2) * sizeof(Q));
            r.c[0] = qi(0);
            for (int i = 0; i <= p.deg; i++) r.c[i + 1] = q_div(p.c[i], qi(i + 1));
        }
        if (r.deg < 0) return vq(qi(0));
        return vpoly(r);
    }
    case N_ROOT: {
        if (!n->c) longjmp(need_series, 1);
        Val l = eval(n->a), r = eval(n->b), g = eval(n->c);
        if (!is_exact(l) || !is_exact(r)) nm_fail("the equation must have exact coefficients");
        if (g.kind != V_Q) nm_fail("the starting value after 'near' must be an exact number");
        Poly p = p_sub(as_poly(l), as_poly(r));
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

/* A quantity with its moment: a + b·o, where o·o is rejected (Methodus, Problem 1). The moment carries the
 * derivative with respect to the unknown term, which is what resolution needs. */
typedef struct { Ser a, b; int hasb; } Dual;

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
        Ser xb = x.hasb ? x.b : s_const(qi(0), x.a.n), yb = y.hasb ? y.b : s_const(qi(0), y.a.n);
        d.b = sign > 0 ? s_add(xb, yb) : s_sub(xb, yb);
    }
    return d;
}

static Dual dmul(Dual x, Dual y) {
    Dual d;
    d.a = s_mul(x.a, y.a);
    d.hasb = x.hasb || y.hasb;
    if (d.hasb) {
        Ser t = s_const(qi(0), d.a.n);
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
        Ser t = x.hasb ? x.b : s_const(qi(0), x.a.n);
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
        d.b = s_mul(s_scale(dfa, alpha), x.b);
    }
    return d;
}

static Dual sev(Node *n, SCtx *cx);
static Ser resolve(Node *root, const char *var, int n);

static Q const_of(Dual d, const char *what) {
    for (int i = 1; i < d.a.n; i++)
        if (q_sign(d.a.c[i])) nm_fail("%s must be a number, not a series", what);
    if (d.hasb) nm_fail("%s must not involve the unknown", what);
    return d.a.n ? d.a.c[0] : qi(0);
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
    if (same(s, len, cx->var, strlen(cx->var))) return dconst(s_var(cx->n));
    Binding *b = lookup(s, len);
    if (!b) nm_fail("the letter %.*s has no value (one letter per series in this version)", (int)len, s);
    switch (b->v.kind) {
    case V_Q: return dconst(s_const(b->v.q, cx->n));
    case V_POLY:
        if (strcmp(b->v.p.var, cx->var)) nm_fail("%s is a polynomial in %s, not in %s", b->name, b->v.p.var, cx->var);
        return dconst(s_from_poly(b->v.p, cx->n));
    case V_REC: return dconst(recipe_series(b, cx->var, cx->n));
    default: nm_fail("%s is an irrational number; irrational coefficients come later", b->name);
    }
    return dconst(s_const(qi(0), cx->n));
}

static Dual sev(Node *n, SCtx *cx) {
    switch (n->kind) {
    case N_NUM: return dconst(s_const(parse_number(n->s, n->len), cx->n));
    case N_NAME: return sname(n->s, n->len, n->primes, cx);
    case N_NEG: return dadd(dconst(s_const(qi(0), cx->n)), sev(n->a, cx), -1);
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
            if (n->kind == N_DERIV) return dconst(s_const(qi(0), cx->n));
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
        if (n->c) nm_fail("a number found with 'near' cannot be a series coefficient unless it is rational");
        return dconst(resolve(n, cx->var, cx->n));
    }
    nm_fail("internal: unknown node");
    return dconst(s_const(qi(0), 0));
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
    Q *c = arena_alloc((size_t)(n + r + 1) * sizeof(Q));
    for (int i = 0; i < n + r + 1; i++) c[i] = qi(0);
    int *given = arena_alloc((size_t)need * sizeof(int));
    memset(given, 0, (size_t)need * sizeof(int));
    SCtx k0; memset(&k0, 0, sizeof k0); k0.var = var; k0.n = 1;
    for (int i = 0; i < root->nconds; i++) {
        if (q_sign(const_of(sev(cs[i].at, &k0), "the point of a start")) != 0)
            nm_fail("starts must be given at 0 in this version");
        int j = cs[i].primes;
        if (j >= need) nm_fail("the equation has order %d; %.*s with %d marks is not a start for it", r, (int)ul, u, j);
        Q v = const_of(sev(cs[i].val, &k0), "a starting value");
        Q fact = qi(1);
        for (int f = 2; f <= j; f++) fact = q_mul(fact, qi(f));
        c[j] = q_div(v, fact);
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
        Ser Y; Y.n = m; Y.c = arena_alloc((size_t)m * sizeof(Q));
        for (int i = 0; i < m; i++) Y.c[i] = i < d ? c[i] : qi(0);
        Ser O = s_const(qi(0), m);
        O.c[d] = qi(1);
        for (int j = 0; j <= r; j++) {
            cx.yd[j].a = s_trunc(Y, k + 1); cx.yd[j].b = s_trunc(O, k + 1); cx.yd[j].hasb = 1;
            Y = s_deriv(Y); O = s_deriv(O);
        }
        cx.n = k + 1;
        Dual e = sev(eqn, &cx);
        Q F = e.a.n > k ? e.a.c[k] : qi(0);
        Q G = (e.hasb && e.b.n > k) ? e.b.c[k] : qi(0);
        if (q_sign(G) == 0)
            nm_fail("the start is a multiple root: the next term is not fixed by the lowest terms (Newton's parallelogram comes later)");
        c[d] = q_neg(q_div(F, G));
    }
    /* substitute back */
    Ser Y; Y.n = n; Y.c = c;
    int check = n - r;
    for (int j = 0; j <= r; j++) { cx.yd[j] = dconst(s_trunc(Y, check)); Y = s_deriv(Y); }
    cx.n = check;
    Dual e = sev(eqn, &cx);
    for (int i = 0; i < e.a.n && i < check; i++)
        if (q_sign(e.a.c[i])) {
            if (i == 0 && r == 0) nm_fail("%.*s(0) = %s does not satisfy the equation at %s = 0", (int)ul, u, q_to_str(c[0]), var);
            nm_fail("substituting the series back leaves a remainder at %s^%d: the equation cannot be resolved term by term here", var, i);
        }
    Ser out; out.n = n; out.c = c;
    return out;
}

/* the letter of a series statement: the one unknown letter that is not the unknown of an equation */
static void letters(Node *n, const char **found, size_t *flen, const char *skip, size_t skiplen) {
    if (!n) return;
    const char *s = NULL; size_t len = 0;
    if (n->kind == N_NAME || n->kind == N_APPLY) {
        Binding *b = lookup(n->s, n->len);
        if (b && b->v.kind == V_REC && n->kind == N_NAME) { s = b->v.var; len = strlen(s); }
        else if (b && b->v.kind == V_POLY) { s = b->v.p.var; len = strlen(s); }
        else if (!b && !(skip && same(n->s, n->len, skip, skiplen))) { s = n->s; len = n->len; }
    }
    if ((n->kind == N_DERIV || n->kind == N_INTEG) && n->len) { s = n->s; len = n->len; }
    if (s) {
        if (*found && !same(*found, *flen, s, len))
            nm_fail("two letters, %.*s and %.*s: a series has one letter in this version", (int)*flen, *found, (int)len, s);
        *found = s; *flen = len;
    }
    if (n->kind == N_ROOT && n->nconds) { skip = n->conds[0].s; skiplen = n->conds[0].len; }
    letters(n->a, found, flen, skip, skiplen);
    letters(n->b, found, flen, skip, skiplen);
    letters(n->c, found, flen, skip, skiplen);
}

/* ---------------- the value of a series at a number ---------------- */

/* The term rule. If the equation is  sum_j p_j(x) y^(j) + h(x) = 0  with polynomials p_j and h, then putting
 * y = sum c_k x^k gives, for n large enough,
 *     D(n) c_n = sum_{t>=1} N_t(n) c_{n-t} - h_{n-s},
 * each term from the ones before, as Newton's A, B, C, D. */
typedef struct { int s, T; Poly D; Poly *N; Poly h; } Rule;

static Poly poly_of(Val v, const char *var) {
    if (v.kind == V_Q) return p_const(v.q);
    if (v.kind == V_POLY && !strcmp(v.p.var, var)) return v.p;
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

static Val series_value(Binding *b, Q a) {
    Node *root = parse_text(b->v.src);
    if (root->kind != N_ROOT || !root->nconds)
        nm_fail("the value of %s needs its definition as an equation with a start (root of ..., y(0) = ...)", b->name);
    const char *u = root->conds[0].s; size_t ul = root->conds[0].len;
    int r = max_primes(root->a, u, ul), r2 = max_primes(root->b, u, ul);
    if (r2 > r) r = r2;
    Rule R = term_rule(root, b->v.var, r);
    if (R.D.deg < 0) nm_fail("the equation fixes no term of %s", b->name);
    Q rr = qabs(a);
    /* how fast the terms shrink in the end: sum over t of lim |N_t / D| r^t must be below 1 */
    Q ginf = qi(0);
    for (int t = 1; t <= R.T; t++) {
        if (R.N[t].deg > R.D.deg) nm_fail("the terms of %s grow too fast: the series has no value away from 0", b->name);
        if (R.N[t].deg == R.D.deg) ginf = q_add(ginf, q_mul(qabs(q_div(R.N[t].c[R.D.deg], R.D.c[R.D.deg])), q_pow(rr, t)));
    }
    if (q_cmp_one(ginf) >= 0)
        nm_fail("at %s the terms of %s shrink too slowly to bound the rest; use smaller arguments and exact relations, as Newton did with 0.1 and 0.2",
                q_to_str(a), b->name);
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
    for (int i = 0; i < nseed; i++) c[i] = seed.c[i];
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
    if (q_sign(a) == 0) return vq(c[0]);
    /* sum term by term in places (the coefficients exact, the powers of a carried as balls); stop when the rest
     * is certainly below the places asked */
    int64_t want = work_prec - GUARD_DIGITS + 10;
    Q eps = q_make(z_from_i64(1), z_pow10(want));
    Ball S = b_from_q(qi(0), work_prec), A = b_from_q(a, work_prec), P = b_from_q(qi(1), work_prec);
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

/* ---------------- writing results ---------------- */

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
            else snprintf(out, 4096, "[exact value %s, rounded to %lld places]", q_to_str(v.q), (long long)places);
        }
        break;
    case V_POLY:
        body = p_to_str(v.p);
        snprintf(out, 4096, "[exact]");
        break;
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
    case V_POLY: {
        Poly p = v.p;
        Q *c = perm_alloc((size_t)(p.deg + 1) * sizeof(Q));
        for (int i = 0; i <= p.deg; i++) c[i] = q_persist(p.c[i]);
        v.p.c = c; v.p.var = dup_perm(p.var, strlen(p.var));
        break;
    }
    case V_BALL: v.ball.m = z_persist(v.ball.m); v.ball.r = z_persist(v.ball.r); break;
    default: break;                                     /* roots and recipes already live in permanent memory */
    }
    return v;
}

static void bind(const char *name, size_t len, Val v) {
    Binding *b = perm_alloc(sizeof *b);
    b->name = dup_perm(name, len);
    b->v = persist_val(v);
    b->next = names; names = b;
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
        expect_op('=');
    }
    const char *src_start = peek()->at;
    Node *e = expr();
    const char *src_end = peek()->at;
    int64_t places = DEFAULT_PLACES, order = DEFAULT_ORDER;
    volatile int asked = 0;
    if (at_key("to")) {
        pos++;
        if (peek()->kind == T_NAME) {                                  /* to x^8 */
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
    if (!setjmp(need_series)) {
        v = eval(e);
        text = show(v, places, asked);
        if (v.kind == V_BALL) {                        /* keep the recipe, so a later request can carry it further */
            memset(&v, 0, sizeof v);
            v.kind = V_NUMREC;
            v.src = dup_perm(src_start, (size_t)(src_end - src_start));
        }
    } else {
        const char *var = NULL; size_t vlen = 0;
        letters(e, &var, &vlen, NULL, 0);
        if (!var) { var = "x"; vlen = 1; }
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
    if (let_name) {
        bind(let_name, let_len, v);
        char *o = arena_alloc(strlen(text) + let_len + 4);
        sprintf(o, "%.*s = %s", (int)let_len, let_name, text);
        return o;
    }
    return text;
}
