CC ?= cc
CFLAGS ?= -O2 -std=c11 -D_POSIX_C_SOURCE=200809L -Wall -Wextra -Wno-unused-parameter -DNM_LIBDIR=\"$(CURDIR)/lib\"
SRC = src/z.c src/q.c src/coef.c src/series.c src/newton.c src/approx.c src/lang.c

newtonmath: $(SRC) src/main.c src/nm.h
	$(CC) $(CFLAGS) -o $@ $(SRC) src/main.c

test: newtonmath tests/check_z
	sh tests/run.sh
	./tests/check_z

tests/check_z: $(SRC) tests/check_z.c src/nm.h
	$(CC) $(CFLAGS) -Isrc -o $@ $(SRC) tests/check_z.c

clean:
	rm -f newtonmath tests/check_z
