# Viva sheet: CRISPR Guide RNA Designer with a Secure Genome Vault

## The problem in plain words

CRISPR-Cas9 edits DNA where a 20-letter guide matches the DNA right next to a short signal, the PAM (NGG for the common Cas9). A good guide matches its target and nowhere else. Cas9 forgives a few wrong letters, so a bad guide cuts other places too. The computing problem: find every candidate guide on both strands, search the whole genome for near-copies fast, rank guides by safety, and protect the genome file.

## The pipeline in five steps

1. **Load** a FASTA genome (or generate a synthetic one with a planted bad guide).
2. **Find** every 20 letters + NGG site on both strands with KMP, Boyer-Moore and Rabin-Karp, which must agree.
3. **Index** every site by pieces of its 20 letters (pigeonhole seeding) so near-copies can be looked up instead of scanned for.
4. **Score and rank** each guide from its near-copies: heap for top K, quicksort for the order, median of medians for the median.
5. **Protect** the data with the vault (Huffman, RSA wrapping a session seed, keystream), and show how a small key falls to Pollard's rho. CpG islands, LCS and alignment are supporting analyses.

## Algorithms, reasons, complexity

| Module | Algorithm | Why it is here | Complexity |
|---|---|---|---|
| Target finder | KMP | Text pointer never backs up | O(n + m) |
| | Boyer-Moore (bad character) | Jumps over text it can rule out | best O(n/m), worst O(nm) |
| | Rabin-Karp | Rolling hash, verify on hit | O(n + m) expected |
| Off-target | Brute force | Baseline | O(n * 20) per guide |
| | Pigeonhole seed index, CSR hash tables, XOR + popcount | Look up instead of scan | build O(S(k+1)), query reads a few buckets |
| Ranking | Min-heap of size K | Keep the best K without sorting everything | O(N log K) |
| | Quicksort | Order the kept guides | O(N log N) average |
| | Median of medians | Median with a worst-case guarantee | O(N) worst case |
| CpG | Prefix sums | Each 200-letter window in O(1) | O(n) |
| | LCS, two rows | Compare two samples | O(mn) time, O(min(m,n)) space |
| Alignment | Needleman-Wunsch | Optimal global score and alignment | O(mn) time and space |
| | Hirschberg | Same score in linear space | O(mn) time, O(n) space |
| | Smith-Waterman (score) | Best local match | O(mn) time, O(n) space |
| Vault | Huffman | Compress before encrypting | O(n + s log s) |
| | Miller-Rabin, modular exponentiation, extended Euclid | Build an RSA key pair | O(log e) multiplications per exponentiation |
| | Pollard's rho | Break a toy key by factoring | about 2^(bits/2) steps |

## Measured results (4.6 Mb synthetic genome unless stated)

- 308,404 targets, identical from all three matchers. On a 16-letter pattern: Boyer-Moore 6.8 ms, KMP 18.6 ms, Rabin-Karp 39.5 ms.
- Off-target, k = 3: brute force 27.4 ms per guide, index 0.080 ms (341x). Index build 54 ms once. With the build counted over 20 queries: 9.8x. 20 of 20 guides gave identical hit sets.
- Planted guide: 1 intended site plus 12 decoys. It ranks worst in its region (specificity 28.16).
- Alignment 8000 x 8000: both methods score 3586. Full table 128 MB, Hirschberg 96 KB.
- Pollard's rho: 220 of 220 keys recovered. Steps grew 1,117x from 12 to 32 bit primes (square-root growth predicts 1,024x).

## Ten likely questions

1. **Why does pigeonhole seeding find every near-match?** Split the guide into k+1 pieces. A site with at most k mismatches has k wrong letters, which can touch at most k pieces. So at least one piece is untouched and matches exactly, and the index finds the site through that piece. Nothing is missed.

2. **How do you avoid reporting the same site twice?** A site can match exactly in several pieces. It is reported only from the first piece that matches exactly. That needs a check against earlier pieces, no hash set.

3. **Why Hirschberg when Needleman-Wunsch already works?** The full table is m times n cells. At 8000 x 8000 that is 128 MB, against 96 KB for Hirschberg. It runs the table fill forward over the first half and backward over the second half, splits where the scores add up best, and recurses. The cost is roughly twice the time.

4. **Why median of medians?** Quickselect can degrade to O(N squared) with bad pivots. Taking the median of group-of-5 medians guarantees a pivot that discards at least about 30% of the items each round, so the total is O(N) in the worst case.

5. **Why a heap for the top K instead of sorting everything?** A size-K heap costs O(N log K) and only the K survivors are sorted at the end. Sorting all N costs O(N log N). With thousands of candidates and K = 10 it saves work and memory.

6. **Why is the vault not secure?** The modulus is at most 64 bits, so rho factors it in milliseconds. The session seed is only 32 bits, so the keystream can be guessed by trying all seeds. The cipher is a simple xorshift generator, not a vetted one. It teaches the pipeline and the attack, nothing more.

7. **Why does Pollard's rho need about 2^(b/2) steps?** It walks a pseudo-random sequence modulo n. Modulo the smaller prime p the sequence repeats after about the square root of p steps (the birthday paradox), and a gcd then reveals p. Square root of 2^b is 2^(b/2), so two more bits of prime means twice the work.

8. **Why does Boyer-Moore win on the 16-letter pattern but not on the 2-letter one?** Its jumps are at most the pattern length. With a 2-letter pattern there is nothing to jump over. With 16 letters it often skips most of the text, which is the O(n/m) best case.

9. **Why are prefix sums used for CpG islands?** Each 200-letter window needs counts of C, G and CG. With prefix sums every count is two subtractions, so sliding the window across the genome is O(n) instead of O(n times 200).

10. **How do you know the index gives the same answer as brute force?** Both return the full list of hits. The program compares the counts, an order-independent signature, and for the tests every hit's position, strand, mismatch count and mismatch mask, on 20 guides and for k from 0 to 5. The planted guide must return exactly 13 sites at k = 3.
