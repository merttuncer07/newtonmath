/* newtonmath: run statements from a file, from -e, or line by line at the terminal. */
#include "nm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

char *nm_run(const char *line, int *failed);

static int run_stream(FILE *in, int interactive) {
    char *line = NULL;
    size_t cap = 0;
    int errors = 0;
    for (;;) {
        if (interactive) { fputs("> ", stdout); fflush(stdout); }
        if (getline(&line, &cap, in) < 0) break;
        int failed;
        char *out = nm_run(line, &failed);
        if (out) puts(out);
        errors += failed;
        arena_reset();
    }
    free(line);
    if (interactive) putchar('\n');
    return errors;
}

int main(int argc, char **argv) {
    if (argc == 3 && !strcmp(argv[1], "-e")) {
        int failed;
        char *out = nm_run(argv[2], &failed);
        if (out) puts(out);
        return failed;
    }
    if (argc == 2) {
        FILE *f = fopen(argv[1], "r");
        if (!f) { perror(argv[1]); return 2; }
        int errors = run_stream(f, 0);
        fclose(f);
        return errors ? 1 : 0;
    }
    if (argc == 1) return run_stream(stdin, isatty(0)) ? 1 : 0;
    fprintf(stderr, "usage: newtonmath [file.nm | -e \"statement\"]\n");
    return 2;
}
