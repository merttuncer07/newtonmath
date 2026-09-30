/* The language: reading a statement, working it out, and writing the result with its verdict.
 *
 *   statement := "let" NAME "=" expr [to]  |  expr [to]
 *   to        := "to" INTEGER "places"
 *   expr      := term { ("+" | "-") term }
 *   term      := unary { ("*" | "/") unary | unary }          (writing side by side multiplies: 2y, 3(x+1))
 *   unary     := "-" unary | power
 *   power     := primary [ "^" unary ]
 *   primary   := NUMBER | NAME | "(" expr ")" | "sqrt" "(" expr ")" | "root" "of" expr "=" expr "near" expr
 *
 * Plain ASCII is the canonical form; the reader also accepts √ × · − ² ³. */
#include "nm.h"

#include <ctype.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern jmp_buf nm_on_error;
extern char nm_error_msg[512];

#define DEFAULT_PLACES 20
#define GUARD_DIGITS 30

/* ---------------- tokens ---------------- */

enum { T_END, T_NUM, T_NAME, T_OP, T_KEY };
typedef struct { int kind; const char *s; size_t len; char op; } Tok;

static const char *KEYWORDS[] = {"let", "to", "places", "root", "of", "near", "sqrt", NULL};

static Tok *toks;
static int ntok, pos;

static int is_key(const char *s, size_t len) {
    for (int i = 0; KEYWORDS[i]; i++)
        if (strlen(KEYWORDS[i]) == len && strncmp(KEYWORDS[i], s, len) == 0) return 1;
    return 0;
}

static void push(int kind, const char *s, size_t len, char op) {
    Tok t = {kind, s, len, op};
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
            i = j; continue;
        }
        if (strchr("+-*/^()=", c)) { push(T_OP, src + i, 1, (char)c); i++; continue; }
        /* keyboard-free aliases (UTF-8) */
        if (!strncmp(src + i, "\xE2\x88\x9A", 3)) { push(T_KEY, "sqrt", 4, 0); i += 3; continue; }      /* √ */
        if (!strncmp(src + i, "\xC3\x97", 2) || !strncmp(src + i, "\xC2\xB7", 2)) {                     /* × · */
            push(T_OP, "*", 1, '*'); i += 2; continue;
        }
        if (!strncmp(src + i, "\xE2\x88\x92", 3)) { push(T_OP, "-", 1, '-'); i += 3; continue; }         /* − */
        if (!strncmp(src + i, "\xC2\xB2", 2) || !strncmp(src + i, "\xC2\xB3", 2)) {                     /* ² ³ */
            push(T_OP, "^", 1, '^');
            push(T_NUM, src[i + 1] == '\xB2' ? "2" : "3", 1, 0);
            i += 2; continue;
        }
        nm_fail("unexpected character '%c'", c);
    }
    push(T_END, "", 0, 0);
}

static Tok *peek(void) { return &toks[pos]; }
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

enum { N_NUM, N_NAME, N_NEG, N_BIN, N_SQRT, N_ROOT };
typedef struct Node { int kind; char op; const char *s; size_t len; struct Node *a, *b, *c; } Node;

static Node *mk(int kind) { Node *n = arena_alloc(sizeof *n); memset(n, 0, sizeof *n); n->kind = kind; return n; }
static Node *bin(char op, Node *a, Node *b) { Node *n = mk(N_BIN); n->op = op; n->a = a; n->b = b; return n; }

static Node *expr(void);
static Node *unary(void);

static Node *primary(void) {
    Tok *t = peek();
    if (t->kind == T_NUM) { Node *n = mk(N_NUM); n->s = t->s; n->len = t->len; pos++; return n; }
    if (t->kind == T_NAME) { Node *n = mk(N_NAME); n->s = t->s; n->len = t->len; pos++; return n; }
    if (at_op('(')) { pos++; Node *e = expr(); expect_op(')'); return e; }
    if (at_key("sqrt")) {
        pos++;
        Node *n = mk(N_SQRT);
        if (at_op('(')) { pos++; n->a = expr(); expect_op(')'); }
        else n->a = primary();                          /* √2 */
        return n;
    }
    if (at_key("root")) {
        pos++; expect_key("of");
        Node *n = mk(N_ROOT);
        n->a = expr(); expect_op('='); n->b = expr(); expect_key("near"); n->c = expr();
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
        else if (peek()->kind == T_NAME || peek()->kind == T_NUM || at_op('(') || at_key("sqrt") || at_key("root"))
            a = bin('*', a, power());                   /* side by side: 2y, 3(x+1) */
        else return a;
    }
}

static Node *expr(void) {
    Node *a = term();
    while (at_op('+') || at_op('-')) { char op = peek()->op; pos++; a = bin(op, a, term()); }
    return a;
}

/* ---------------- values ---------------- */

enum { V_Q, V_POLY, V_ROOT, V_BALL };
typedef struct { int kind; Q q; Poly p; Root *root; Ball ball; } Val;

typedef struct Binding { char *name; Val v; struct Binding *next; } Binding;
static Binding *names;

static Val vq(Q q) { Val v; memset(&v, 0, sizeof v); v.kind = V_Q; v.q = q; return v; }
static Val vpoly(Poly p) {
    if (p.deg <= 0) return vq(p.deg < 0 ? q_from_z(z_zero()) : p.c[0]);
    Val v; memset(&v, 0, sizeof v); v.kind = V_POLY; v.p = p; return v;
}
static Val vball(Ball b) { Val v; memset(&v, 0, sizeof v); v.kind = V_BALL; v.ball = b; return v; }

static int64_t work_prec;                               /* significant digits for balls in this statement */

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

static Q parse_number(const char *s, size_t len) {
    const char *dot = memchr(s, '.', len);
    if (!dot) return q_from_z(z_from_dec(s, len));
    size_t ip = (size_t)(dot - s), fp = len - ip - 1;
    char *all = arena_alloc(len + 2);
    memcpy(all, s, ip); memcpy(all + ip, dot + 1, fp);
    if (ip + fp == 0) { all[0] = '0'; fp = 0; ip = 1; }
    return q_make(z_from_dec(all, ip + fp), z_pow10((int64_t)fp));   /* 0.1 is exactly 1/10 */
}

/* c^(1/n) for a rational c: exact when c is a perfect n-th power, otherwise a root of B*y^n - A */
static Val nth_root(Q c, int64_t n) {
    if (n < 2 || n > 1000) nm_fail("root index out of range");
    int neg = q_sign(c) < 0;
    if (neg && n % 2 == 0) nm_fail("an even root of a negative number is not real (complex numbers come later)");
    Z A = z_abs(c.num), B = c.den;
    Z ra = z_iroot(A, (unsigned)n), rb = z_iroot(B, (unsigned)n);
    if (z_cmp(z_pow(ra, (unsigned)n), A) == 0 && z_cmp(z_pow(rb, (unsigned)n), B) == 0) {
        Q r = q_make(ra, rb);
        return vq(neg ? q_neg(r) : r);
    }
    /* starting value from the integer root of the scaled number, as in root extraction by hand */
    int64_t D0 = 6;
    Z scaled = z_iroot(z_div_round(z_mul_pow10(A, D0 * n), B), (unsigned)n);
    Q guess = q_make(neg ? z_neg(scaled) : scaled, z_pow10(D0));
    Poly p = p_var("y");
    p = p_pow(p, (unsigned)n);
    p = p_sub(p_scale(p, q_from_z(B)), p_const(q_from_z(neg ? z_neg(A) : A)));
    Val v; memset(&v, 0, sizeof v);
    v.kind = V_ROOT;
    v.root = root_new(p, guess, D0);
    return v;
}

static Val eval(Node *n);

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
        if (den != 1 || num < 0) nm_fail("a letter can only be raised to a whole non-negative power");
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
            if (b.kind != V_Q) nm_fail("dividing by a letter is not in this version");
            return vpoly(p_scale(x, q_div(q_from_z(z_from_i64(1)), b.q)));
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

static Val eval(Node *n) {
    switch (n->kind) {
    case N_NUM: return vq(parse_number(n->s, n->len));
    case N_NAME: {
        for (Binding *b = names; b; b = b->next)
            if (strlen(b->name) == n->len && !strncmp(b->name, n->s, n->len)) return b->v;
        char *name = arena_alloc(n->len + 1);
        memcpy(name, n->s, n->len); name[n->len] = 0;
        return vpoly(p_var(name));                      /* an unknown letter */
    }
    case N_NEG: {
        Val v = eval(n->a);
        return arith('-', vq(q_from_z(z_zero())), v);
    }
    case N_BIN: return arith(n->op, eval(n->a), eval(n->b));
    case N_SQRT: return power_val(eval(n->a), vq(q_make(z_from_i64(1), z_from_i64(2))));
    case N_ROOT: {
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
    return vq(q_from_z(z_zero()));
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
        if (!r->certified) {
            body = fixed_str(z_div_round(r->X, z_pow10(r->D - places)), places);
            snprintf(out, 4096, "[not certified: %s shows no sign change here; a double root?]", root_equation_str(r));
            break;
        }
        body = fixed_str(z_div_round(r->X, z_pow10(r->D - places)), places);
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
        size_t l = strlen(p.var);
        char *name = perm_alloc(l + 1);
        memcpy(name, p.var, l + 1);
        v.p.c = c; v.p.var = name;
        break;
    }
    case V_BALL: v.ball.m = z_persist(v.ball.m); v.ball.r = z_persist(v.ball.r); break;
    case V_ROOT: break;                                  /* roots already live in permanent memory */
    }
    return v;
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
    const char *let_name = NULL; size_t let_len = 0;
    if (at_key("let")) {
        pos++;
        if (peek()->kind != T_NAME) nm_fail("expected a name after 'let'");
        let_name = peek()->s; let_len = peek()->len; pos++;
        expect_op('=');
    }
    Node *e = expr();
    int64_t places = DEFAULT_PLACES;
    int asked = 0;
    if (at_key("to")) {
        pos++;
        if (peek()->kind != T_NUM || memchr(peek()->s, '.', peek()->len)) nm_fail("expected a whole number after 'to'");
        Z k = z_from_dec(peek()->s, peek()->len);
        if (!z_fits_i64(k, &places) || places > 100000) nm_fail("too many places");
        pos++;
        expect_key("places");
        asked = 1;
    }
    if (peek()->kind != T_END) nm_fail("unexpected '%.*s'", (int)peek()->len, peek()->s);
    work_prec = places + GUARD_DIGITS;
    Val v = eval(e);
    char *text = show(v, places, asked);
    if (let_name) {
        Binding *b = perm_alloc(sizeof *b);
        b->name = perm_alloc(let_len + 1);
        memcpy(b->name, let_name, let_len); b->name[let_len] = 0;
        b->v = persist_val(v);
        b->next = names; names = b;
        char *o = arena_alloc(strlen(text) + let_len + 4);
        sprintf(o, "%s = %s", b->name, text);
        return o;
    }
    return text;
}
