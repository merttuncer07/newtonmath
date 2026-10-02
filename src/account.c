/* The account a statement gives of itself: the steps actually taken ("show") and a few facts the computation
 * produced anyway. Written to the terminal for a person, or kept as JSON for a program (newtonmath -j). Nothing is
 * computed here; the engine only reports what it did. */
#include "nm.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int nm_show_work, nm_work_lines, nm_json;

typedef struct { char *s; size_t n, cap; } Buf;
static Buf work, facts;

static void put(Buf *b, const char *s, size_t n) {
    if (b->n + n + 1 > b->cap) {
        b->cap = (b->n + n + 1) * 2;
        b->s = realloc(b->s, b->cap);
        if (!b->s) { fprintf(stderr, "out of memory\n"); exit(2); }
    }
    memcpy(b->s + b->n, s, n);
    b->n += n;
    b->s[b->n] = 0;
}

static void put_json_str(Buf *b, const char *s) {
    put(b, "\"", 1);
    for (; *s; s++) {
        char e[8];
        if (*s == '"' || *s == '\\') { e[0] = '\\'; e[1] = *s; put(b, e, 2); }
        else if (*s == '\n') put(b, "\\n", 2);
        else if ((unsigned char)*s < 0x20) { snprintf(e, sizeof e, "\\u%04x", *s); put(b, e, 6); }
        else put(b, s, 1);
    }
    put(b, "\"", 1);
}

static char *vfmt(const char *fmt, va_list ap) {
    va_list cp;
    va_copy(cp, ap);
    int n = vsnprintf(NULL, 0, fmt, cp);
    va_end(cp);
    char *s = arena_alloc((size_t)n + 1);
    vsnprintf(s, (size_t)n + 1, fmt, ap);
    return s;
}

void nm_account_reset(void) { work.n = facts.n = 0; if (work.s) work.s[0] = 0; if (facts.s) facts.s[0] = 0; nm_work_lines = 0; }

void nm_work(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    char *line = vfmt(fmt, ap);
    va_end(ap);
    if (nm_json) { if (work.n) put(&work, ",", 1); put_json_str(&work, line); }
    else printf("  %s\n", line);
    nm_work_lines++;
}

/* a fact: `json` is already a JSON value (a number, true/false, an array, or a quoted string) */
void nm_fact(const char *key, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    char *val = vfmt(fmt, ap);
    va_end(ap);
    Buf k = {0, 0, 0};                           /* a key is given once per statement: the first one stands */
    put_json_str(&k, key);
    put(&k, ":", 1);
    int dup = facts.n && strstr(facts.s, k.s) != NULL;
    free(k.s);
    if (dup) return;
    if (facts.n) put(&facts, ",", 1);
    put_json_str(&facts, key);
    put(&facts, ":", 1);
    put(&facts, val, strlen(val));
}

void nm_account_print_facts(void) { if (facts.n && !nm_json) printf("  facts: {%s}\n", facts.s); }

char *nm_json_str(const char *s) {               /* s as a JSON string, in the arena */
    Buf b = {0, 0, 0};
    put_json_str(&b, s);
    char *r = arena_alloc(b.n + 1);
    memcpy(r, b.s, b.n + 1);
    free(b.s);
    return r;
}

/* how sure the answer is, from its verdict: the first word, as the verdicts are written */
static const char *status_of(const char *verdict) {
    static const char *words[] = {"proved", "certified", "bounded", "exact", "probable", "rule", NULL};
    if (strstr(verdict, "probable")) return "probable";
    if (strstr(verdict, "proved")) return "proved";
    for (int i = 0; words[i]; i++) if (!strncmp(verdict, words[i], strlen(words[i]))) return words[i];
    if (strstr(verdict, "certified")) return "certified";
    if (strstr(verdict, ": exact")) return "exact";
    return "unlabelled";
}

/* one line of JSON for a statement and its result */
char *nm_account_json(const char *input, const char *out, int failed) {
    Buf b = {0, 0, 0};
    put(&b, "{\"input\":", 9);
    char *in = arena_alloc(strlen(input) + 1);
    strcpy(in, input);
    size_t n = strlen(in);
    while (n && (in[n - 1] == '\n' || in[n - 1] == '\r')) in[--n] = 0;
    put_json_str(&b, in);
    if (failed) {
        put(&b, ",\"error\":", 9);
        put_json_str(&b, !strncmp(out, "error: ", 7) ? out + 7 : out);
    } else {
        const char *last = NULL;                  /* the verdict: the last "  [...]" or "\n[...]" at the end */
        for (const char *p = out; p && (p = strstr(p, "[")); p++)
            if (p - out >= 2 && (p[-1] == '\n' || (p[-1] == ' ' && p[-2] == ' '))) last = p - (p[-1] == '\n' ? 1 : 2);
        size_t alen = out ? strlen(out) : 0;
        if (last && out[alen - 1] == ']') alen = (size_t)(last - out);
        char *ans = arena_alloc(alen + 1);
        if (out) memcpy(ans, out, alen);
        ans[alen] = 0;
        put(&b, ",\"answer\":", 10);
        put_json_str(&b, ans);
        if (alen < (out ? strlen(out) : 0)) {
            size_t skip = out[alen] == '\n' ? 2 : 3, vl = strlen(out) - alen - skip - 1;
            char *verdict = arena_alloc(vl + 1);
            memcpy(verdict, out + alen + skip, vl);
            verdict[vl] = 0;
            put(&b, ",\"status\":", 10);
            put_json_str(&b, status_of(verdict));
            put(&b, ",\"verdict\":", 11);
            put_json_str(&b, verdict);
        } else if (strstr(ans, "choose one with")) {   /* several branches: the caller must pick a start */
            put(&b, ",\"status\":\"needs_choice\"", 24);
        } else if (strstr(ans, "O(")) {           /* a series: every coefficient shown is exact, the rest is O(...) */
            put(&b, ",\"status\":\"exact_to_order\"", 26);
        }
    }
    if (facts.n) { put(&b, ",\"facts\":{", 10); put(&b, facts.s, facts.n); put(&b, "}", 1); }
    if (work.n) { put(&b, ",\"work\":[", 9); put(&b, work.s, work.n); put(&b, "]", 1); }
    put(&b, "}", 1);
    char *r = arena_alloc(b.n + 1);
    memcpy(r, b.s, b.n + 1);
    free(b.s);
    return r;
}
