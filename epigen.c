/* epigen.c: CpG islands and sequence comparison.
 *
 * CpG islands: prefix sums of C, G and CpG counts make every 200-letter window
 * an O(1) lookup (two subtractions per count). A window qualifies when
 * GC > 50% and observed/expected CpG > 0.6, where
 *     expected = (C count * G count) / window length.
 * Overlapping or touching qualifying windows are merged into one island.
 * Whole scan: O(n) time, O(n) memory for the three prefix arrays.
 *
 * LCS (longest common subsequence): dynamic programming over a table, but only
 * two rows are kept. O(m*n) time, O(min(m,n)) space. Gives the length only.
 */
#include "crispr.h"
#include <stdlib.h>
#include <string.h>

void islandvec_free(IslandVec *v) { free(v->a); v->a = NULL; v->n = v->cap = 0; }

static void island_push(IslandVec *v, CpgIsland x) {
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 16;
        v->a = realloc(v->a, v->cap * sizeof *v->a);
        if (!v->a) exit(1);
    }
    v->a[v->n++] = x;
}

void cpg_islands(const char *seq, size_t n, IslandVec *out) {
    const size_t W = 200;
    if (n < W) return;
    uint32_t *pc = xcalloc(n + 1, sizeof *pc);
    uint32_t *pg = xcalloc(n + 1, sizeof *pg);
    uint32_t *pcg = xcalloc(n + 1, sizeof *pcg);
    for (size_t i = 0; i < n; i++) {
        pc[i + 1] = pc[i] + (seq[i] == 'C');
        pg[i + 1] = pg[i] + (seq[i] == 'G');
        pcg[i + 1] = pcg[i] + (i + 1 < n && seq[i] == 'C' && seq[i + 1] == 'G');
    }
    int open = 0;
    size_t cs = 0, ce = 0;
    for (size_t s = 0; s + W <= n; s++) {
        double c = pc[s + W] - pc[s], g = pg[s + W] - pg[s];
        double cg = pcg[s + W - 1] - pcg[s]; /* dinucleotides fully inside the window */
        int pass = (c + g) / W > 0.5 && c > 0 && g > 0 && cg * W / (c * g) > 0.6;
        if (!pass) continue;
        if (open && s <= ce) { ce = s + W; continue; }
        if (open) {
            double cc = pc[ce] - pc[cs], gg = pg[ce] - pg[cs], x = pcg[ce - 1] - pcg[cs];
            island_push(out, (CpgIsland){cs, ce, (cc + gg) / (double)(ce - cs), x * (double)(ce - cs) / (cc * gg)});
        }
        open = 1; cs = s; ce = s + W;
    }
    if (open) {
        double cc = pc[ce] - pc[cs], gg = pg[ce] - pg[cs], x = pcg[ce - 1] - pcg[cs];
        island_push(out, (CpgIsland){cs, ce, (cc + gg) / (double)(ce - cs), x * (double)(ce - cs) / (cc * gg)});
    }
    free(pc); free(pg); free(pcg);
}

void gc_track(const char *seq, size_t n, size_t bins, double *out) {
    for (size_t b = 0; b < bins; b++) {
        size_t lo = b * n / bins, hi = (b + 1) * n / bins, gc = 0;
        for (size_t i = lo; i < hi; i++) gc += (seq[i] == 'G' || seq[i] == 'C');
        out[b] = hi > lo ? (double)gc / (double)(hi - lo) : 0.0;
    }
}

size_t lcs_length(const char *a, size_t la, const char *b, size_t lb) {
    if (lb > la) { const char *t = a; a = b; b = t; size_t u = la; la = lb; lb = u; }
    size_t *prev = xcalloc(lb + 1, sizeof *prev), *cur = xcalloc(lb + 1, sizeof *cur);
    for (size_t i = 1; i <= la; i++) {
        for (size_t j = 1; j <= lb; j++) {
            if (a[i - 1] == b[j - 1]) cur[j] = prev[j - 1] + 1;
            else cur[j] = prev[j] > cur[j - 1] ? prev[j] : cur[j - 1];
        }
        size_t *t = prev; prev = cur; cur = t;
    }
    size_t r = prev[lb];
    free(prev); free(cur);
    return r;
}
