CC ?= cc
CFLAGS ?= -O2 -std=c11 -D_POSIX_C_SOURCE=200809L -Wall -Wextra -Wno-unused-parameter
SRC = src/z.c src/q.c src/coef.c src/series.c src/newton.c src/elim.c src/ratfun.c src/integ.c src/sqrtint.c src/linalg.c src/cplx.c src/arith.c src/approx.c src/account.c src/factor.c src/lang.c src/libs.c
OBJ = $(SRC:.c=.o)
CORE = $(filter-out src/lang.o src/libs.o,$(OBJ))

newtonmath: $(OBJ) src/main.o
	$(CC) $(CFLAGS) -o $@ $(OBJ) src/main.o

# the libraries in lib/ are compiled into the program: one file that runs anywhere
src/libs.c: $(wildcard lib/*.nm) Makefile
	@{ echo '/* generated from the files in lib/ by the Makefile; do not edit */'; \
	  echo 'const char *nm_lib_names[] = {'; for f in lib/*.nm; do b=$$(basename $$f .nm); echo "  \"$$b\","; done; echo '  0};'; \
	  echo 'const char *nm_lib_texts[] = {'; for f in lib/*.nm; do awk '{ gsub(/\\/, "\\\\"); gsub(/"/, "\\\""); printf "  \"%s\\n\"\n", $$0 } END { print "  ," }' $$f; done; echo '  0};'; } > $@

src/%.o: src/%.c src/nm.h
	$(CC) $(CFLAGS) -c -o $@ $<

UNITS = $(patsubst tests/unit/%.c,tests/unit/%_test,$(wildcard tests/unit/*.c))

test: newtonmath tests/check_z $(UNITS)
	sh tests/run.sh
	sh tests/roundtrip.sh
	./tests/check_z
	for u in $(UNITS); do ./$$u > /dev/null || { ./$$u; exit 1; }; ./$$u | tail -1; done

tests/unit/%_test: tests/unit/%.c $(CORE) src/nm.h
	$(CC) $(CFLAGS) -Isrc -o $@ $(CORE) $<

tests/check_z: $(OBJ) tests/check_z.c src/nm.h
	$(CC) $(CFLAGS) -Isrc -o $@ $(OBJ) tests/check_z.c

clean:
	rm -f newtonmath src/libs.c tests/check_z tests/unit/*_test src/*.o
