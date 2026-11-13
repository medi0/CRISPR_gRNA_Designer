CC      = gcc
CFLAGS  = -std=c11 -O2 -Wall -Wextra
SRC     = main.c genome.c strmatch.c offtarget.c scoring.c epigen.c align.c vault.c bench.c server.c

crispr: $(SRC) crispr.h
	$(CC) $(CFLAGS) -o crispr $(SRC) -lm

asan: $(SRC) crispr.h
	$(CC) -std=c11 -g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra -o crispr_asan $(SRC) -lm

test: crispr
	./crispr demo
	$(CC) $(CFLAGS) -I. -o test_phase1 tests/test_phase1.c genome.c strmatch.c offtarget.c scoring.c
	./test_phase1
	python3 tests/api_test.py

clean:
	rm -f crispr crispr_asan test_phase1 crispr.exe results.csv

.PHONY: asan test clean
