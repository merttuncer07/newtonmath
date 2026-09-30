CC ?= cc
CFLAGS ?= -O2 -std=c11 -D_POSIX_C_SOURCE=200809L -Wall -Wextra -Wno-unused-parameter -DNM_LIBDIR=\"$(CURDIR)/lib\"
SRC = src/z.c src/q.c src/coef.c src/series.c src/newton.c src/elim.c src/approx.c src/lang.c

newtonmath: $(SRC) src/main.c src/nm.h
	$(CC) $(CFLAGS) -o $@ $(SRC) src/main.c

UNITS = $(patsubst tests/unit/%.c,tests/unit/%_test,$(wildcard tests/unit/*.c))

test: newtonmath tests/check_z $(UNITS)
	sh tests/run.sh
	./tests/check_z
	for u in $(UNITS); do ./$$u > /dev/null || { ./$$u; exit 1; }; ./$$u | tail -1; done

tests/unit/%_test: tests/unit/%.c $(SRC) src/nm.h
	$(CC) $(CFLAGS) -Isrc -o $@ $(SRC:src/lang.c=) $<

tests/check_z: $(SRC) tests/check_z.c src/nm.h
	$(CC) $(CFLAGS) -Isrc -o $@ $(SRC) tests/check_z.c

clean:
	rm -f newtonmath tests/check_z tests/unit/*_test
