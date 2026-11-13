# CRISPR Guide RNA Designer with a Secure Genome Vault

A mini project for a Design and Analysis of Algorithms lab. It picks CRISPR guide RNAs for a genome, checks how safe each one is, and locks the genome file. Every feature exists to show an algorithm, its complexity, and a measured result. It is not a production bioinformatics tool.

All algorithms are written from scratch in C11 with no third-party libraries. The web page is one HTML file served by the same C program.

## Build

Linux:

```
make
```

or without make:

```
gcc -std=c11 -O2 -Wall -Wextra -o crispr main.c genome.c strmatch.c offtarget.c scoring.c epigen.c align.c vault.c bench.c server.c -lm
```

Windows (PowerShell, 64-bit MinGW-w64). First check that this prints `x86_64`:

```
gcc -dumpmachine
```

Then build with one command:

```
gcc -std=gnu11 -O2 -Wall -Wextra -o crispr.exe main.c genome.c strmatch.c offtarget.c scoring.c epigen.c align.c vault.c bench.c server.c -lpsapi -lws2_32
```

`unsigned __int128` needs 64-bit MinGW. **The Windows build has not been run by the author** (see Limitations).

## Run

The fastest way to see everything:

```
./crispr demo
```

Start the web interface, then open the address it prints:

```
./crispr serve 8080
```

Run it from this folder, because it reads `web/index.html` relative to the current directory. On Windows use `.\crispr.exe` in place of `./crispr`.

## Commands

A `<genome>` is a FASTA file, or `synth:N` for a built-in synthetic genome of N bases (N at least 20000).

| Command | What it does |
|---|---|
| `demo [n]` | Runs every module on a synthetic genome (default 300000) and prints PASS/FAIL checks |
| `gen <n> <out.fa> [seed]` | Writes a synthetic genome as FASTA |
| `find <genome>` | Finds all NGG targets with KMP, Boyer-Moore and Rabin-Karp and compares them |
| `offtarget <genome> <guide> [k]` | Every site within k mismatches (0 to 5, default 3), index against brute force |
| `rank <genome> <from> <to> [K] [k]` | Best and worst guides in a region |
| `cpg <genome>` | CpG islands |
| `lcs <genomeA> <genomeB> [len]` | Longest common subsequence of the first len letters |
| `align <A> <B>` | A and B are letters or `rand:LENGTH` (up to 8000) |
| `keygen <bits> <prefix>` | Writes `prefix.pub` and `prefix.priv` (40 to 64 bits) |
| `enc <in> <out> <key.pub>` | Compresses and encrypts a file |
| `dec <in> <out> <key.priv>` | Decrypts and restores a file |
| `attack [minbits maxbits]` | Factors toy keys with Pollard's rho (8 to 32 bit primes) |
| `bench [out.csv]` | Timing runs, optionally written to CSV |
| `serve [port]` | Web interface on 127.0.0.1 only (default 8080) |

Positions are 0-based and mark the start of the 23-letter site (20-letter guide plus 3-letter PAM) on the forward strand.

## Files

| File | Contents |
|---|---|
| `crispr.h` | Shared types and prototypes |
| `genome.c` | FASTA parse/load/save, reverse complement, random numbers, synthetic genome, mutation, dynamic arrays, timer |
| `strmatch.c` | KMP, Boyer-Moore, Rabin-Karp, target finder |
| `offtarget.c` | Seed index, query, brute-force baseline |
| `scoring.c` | Risk scoring, heap top-K, quicksort, median of medians |
| `epigen.c` | CpG islands, LCS |
| `align.c` | Needleman-Wunsch, Hirschberg, Smith-Waterman |
| `vault.c` | Huffman, Miller-Rabin, RSA, stream cipher, Pollard's rho |
| `bench.c` | Timing runs used by `bench` and `/api/bench` (an addition to the planned layout, so both the CLI and the server share one implementation) |
| `server.c` | HTTP server and JSON API |
| `main.c` | Command-line interface |
| `web/index.html` | The whole front end, no CDN, no web fonts, works offline |
| `tests/` | Phase 1 C test, API test (Python), headless DOM test (Node and jsdom) |

Modules return data and do not print. Printing happens in `main.c`, JSON is built only in `server.c`.

## Measured results

Linux x86_64, GCC 13, `-O2`. Synthetic genome, seed 42. Timings vary a few percent from run to run.

**Target finder, 4.6 Mb genome, 308,404 targets.** All three matchers return the same count and the same order-independent signature.

| Pattern | KMP | Boyer-Moore | Rabin-Karp |
|---|---|---|---|
| 2-letter (GG and CC, whole finder) | 34.6 ms | 37.5 ms | 81.3 ms |
| 16-letter, one scan | 18.6 ms | 6.8 ms | 39.5 ms |

**Off-target search, 4.6 Mb, k = 3.**

| Measure | Result |
|---|---|
| Sites indexed | 308,404 (4 segments of 5 letters) |
| Index build | 54.3 ms, once |
| Brute force per query | 27.4 ms |
| Index per query | 0.080 ms (341x faster) |
| 20 queries including the build | 548 ms against 56 ms (9.8x), break-even after 2 queries |
| Identical results | 20 of 20 guides (count, signature and every hit) |
| Planted guide | 13 sites: 1 intended, 4 with 1 mismatch, 4 with 2, 4 with 3 |

The index was also checked against brute force for k = 0 to 5 on four guides each.

**Ranking, 30 kb region around the planted guide on 4.6 Mb.** 2,069 candidates, 207 dropped for TTTT, 1,862 scored in 45.6 ms. The planted guide is the worst one: specificity 28.16 with 12 off-targets. The median specificity is 100, since most guides have no near-copies. The median-of-medians function was also checked against a full sort on 300 random arrays, including ones full of duplicates, with 0 mismatches.

**Alignment, 8000 x 8000 random sequences.** Needleman-Wunsch and Hirschberg both score 3586. The full table takes 128,064,004 bytes and Hirschberg 96,010 bytes (1,334x less). Smith-Waterman local score is 3608 using 32,004 bytes. The process peak was 123 MB, almost all of it the full table. Hirschberg is slower: at length 4000 it took 78.5 ms against 39.4 ms.

**Benchmarks (`./crispr bench`).**

| Genome | KMP | BM | RK | CpG scan | Index build | Index query | Brute query |
|---|---|---|---|---|---|---|---|
| 100,000 | 0.77 | 1.04 | 2.12 | 0.99 | 0.9 | 0.0012 | 0.69 |
| 300,000 | 2.29 | 4.22 | 6.62 | 4.35 | 4.1 | 0.0080 | 2.08 |
| 1,000,000 | 9.11 | 9.50 | 18.86 | 12.55 | 11.1 | 0.0257 | 6.12 |
| 4,600,000 | 40.65 | 40.46 | 85.52 | 54.41 | 59.5 | 0.0918 | 30.95 |

All values in milliseconds. Alignment: length 4000 uses 32,032,004 bytes for the full table and 48,010 bytes for Hirschberg.

**Vault.** A 50,738-byte FASTA file became 14,012 bytes compressed and 14,048 bytes encrypted. The restored file was byte-identical. A different key pair, and the right modulus with a wrong private exponent, were both rejected.

**Pollard's rho.** The private key was recovered in all 220 attempts (20 per size, 11 sizes, primes of 12 to 32 bits). Average steps went from 34 at 12 bits to 38,472 at 32 bits, a growth of 1,117x. Pure 2^(b/2) growth predicts 1,024x. The measured step count sits at roughly 0.6 of the 2^(b/2) line, which is a constant factor and does not change the growth.

## What was tested

| Test | Result |
|---|---|
| Build, Linux, `-Wall -Wextra` | zero warnings |
| Build with `-fsanitize=address,undefined` | zero warnings; `crispr demo` (12 of 12 checks) and the API test run with no sanitizer reports |
| Three matchers agree | pass, 300 kb and 4.6 Mb |
| Index equals brute force | pass, 20 guides at k = 3 on 4.6 Mb and 300 kb |
| Planted guide: 1 intended + 12 decoys at k = 3 | pass, each decoy found at its planted position with its planted mismatch count |
| Planted guide among the worst; guides with no near-matches score 100 | pass |
| Hirschberg score equals full-table score | pass, 8000 x 8000 and in the demo |
| Vault round trip, wrong key rejected | pass, CLI, demo and `/api/vault` |
| Pollard's rho recovers keys, growth about 2^(b/2) | pass |
| Every API endpoint with `curl`-style calls, a FASTA upload, a garbage upload, bad parameters | pass (`tests/api_test.py`) |
| Page loaded in jsdom, every page clicked, error paths checked, `node --check` | pass, with one assertion that jsdom cannot evaluate (below) |

Run `make test` for the C demo, the Phase 1 test and the API test. The DOM test needs `npm install jsdom` in a scratch folder and `node tests/dom_test.js` with the server binary in the current directory.

## Limitations

- **Not tested:** the Windows build; the page in a real browser (only jsdom, which does not compute CSS variables, so the red colour of error text was checked by reading the CSS rule instead); an actual NCBI genome such as E. coli. The FASTA loader handles lower case, N, and Windows line endings and was checked on a short string and on 4.6 Mb synthetic FASTA files, not on a real `.fna` download.
- **The vault is a teaching implementation.** The RSA modulus is at most 64 bits, the session seed is only 32 bits, and the cipher is a simple xorshift keystream. Anyone can break it in seconds. Do not use it for real data.
- **Pollard's rho cannot be run on large keys here.** Infeasibility for real key sizes is shown by the measured growth curve and an extrapolation, not by a failed attack. Rho is also not the best known attack on large keys.
- **The scoring is a teaching version.** The position weights are the published MIT (Hsu 2013) table as the author typed it in; check them against the paper before relying on them. Specificity is 100 / (1 + sum of risks) as specified, which is scaled differently from the MIT web tool.
- Only SpCas9 (NGG) is modelled. Off-targets are substitutions only, with no insertions or deletions in the match.
- The seed index handles k up to 5 only, and the whole index lives in memory.
- Needleman-Wunsch uses 16-bit cells and refuses sequences whose lengths add up to more than 16,000. Web alignment input is capped at 3,000 letters each for pasted text and 8,000 for the random pair.
- The server is single-threaded and blocks while it works, so the benchmark and the 8000-letter alignment hold up other requests for a few seconds. It listens on 127.0.0.1 only and accepts uploads up to 64 MB.
- The synthetic genome is a toy. Real genomes have repeats that make off-target counts much messier.
