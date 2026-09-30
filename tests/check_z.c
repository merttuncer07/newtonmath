/* The integer core checked by inverse operations on random numbers, and against bc as a second route. */
#include "nm.h"

#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern jmp_buf nm_on_error;
extern char nm_error_msg[512];

static unsigned long long seed = 88172645463325252ULL;
static unsigned rnd(void) { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; return (unsigned)seed; }

static Z random_z(int maxdigits) {
    int len = 1 + (int)(rnd() % (unsigned)maxdigits);
    char *s = arena_alloc((size_t)len + 1);
    for (int i = 0; i < len; i++) s[i] = (char)('0' + rnd() % 10);
    /* sometimes long runs of 9s and 0s, where carries and borrows go far */
    if (rnd() % 4 == 0) for (int i = 0; i < len; i++) s[i] = (rnd() % 2) ? '9' : s[i];
    if (rnd() % 4 == 0) for (int i = 1; i < len; i++) s[i] = (rnd() % 2) ? '0' : s[i];
    Z z = z_from_dec(s, (size_t)len);
    return (rnd() % 2) ? z_neg(z) : z;
}

int main(void) {
    int fail = 0, checks = 0;
    if (setjmp(nm_on_error)) { printf("check_z: error: %s\n", nm_error_msg); return 1; }
    for (int t = 0; t < 3000; t++) {
        Z a = random_z(t < 2000 ? 40 : 400), b = random_z(t < 2000 ? 40 : 300);
        if (z_is_zero(b)) continue;
        Z q, r;
        checks += 3;
        if (z_cmp(z_sub(z_add(a, b), b), a)) { printf("add/sub: %s %s\n", z_to_str(a), z_to_str(b)); fail++; }
        z_divmod(z_mul(a, b), b, &q, &r);
        if (z_cmp(q, a) || !z_is_zero(r)) { printf("mul/div: %s %s\n", z_to_str(a), z_to_str(b)); fail++; }
        z_divmod(a, b, &q, &r);                         /* divmod checks q*b + r = a itself */
        Z s = z_iroot(z_abs(a), 2);
        if (z_cmp(z_mul(s, s), z_abs(a)) > 0 || z_cmp(z_mul(z_add(s, z_from_i64(1)), z_add(s, z_from_i64(1))), z_abs(a)) <= 0) {
            printf("isqrt: %s\n", z_to_str(a)); fail++;
        }
        if (t % 50 == 0) arena_reset();
    }
    /* bc as an independent route for products and quotients */
    {
        char path[] = "/tmp/nm_check_XXXXXX";
        int fd = mkstemp(path);
        if (fd >= 0) {
            FILE *in = fdopen(fd, "w");
            Z as[40], bs[40];
            for (int i = 0; i < 40; i++) {
                as[i] = random_z(300); bs[i] = random_z(150);
                if (z_is_zero(bs[i])) bs[i] = z_from_i64(7);
                fprintf(in, "%s*%s\n%s/%s\n", z_to_str(as[i]), z_to_str(bs[i]), z_to_str(as[i]), z_to_str(bs[i]));
            }
            fclose(in);
            char cmd[256];
            snprintf(cmd, sizeof cmd, "BC_LINE_LENGTH=0 bc < %s", path);
            FILE *out = popen(cmd, "r");
            char *line = NULL; size_t cap = 0;
            for (int i = 0; i < 40 && out; i++) {
                Z q, r;
                z_divmod(as[i], bs[i], &q, &r);
                const char *want[2] = {z_to_str(z_mul(as[i], bs[i])), z_to_str(q)};
                for (int k = 0; k < 2; k++) {
                    if (getline(&line, &cap, out) < 0) { printf("bc ended early\n"); fail++; break; }
                    line[strcspn(line, "\n")] = 0;
                    checks++;
                    if (strcmp(line, want[k])) { printf("bc disagrees: %s vs %s\n", line, want[k]); fail++; }
                }
            }
            if (out) pclose(out);
            free(line);
            remove(path);
        }
    }
    printf("check_z: %d checks, %d failed\n", checks, fail);
    return fail != 0;
}
