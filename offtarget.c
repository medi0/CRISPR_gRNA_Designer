/* offtarget.c: find every PAM-adjacent site within k substitutions of a guide.
 *
 * Brute force:  slide the guide over every genome position on both strands and
 *               count mismatches. O(n * L) per guide (L = 20).
 *
 * Seed index:   pigeonhole principle. Split the 20 letters into k+1 segments.
 *               If a site has at most k mismatches, at least one segment has
 *               zero mismatches (k mismatches cannot hit k+1 segments).
 *               So we hash every PAM-adjacent site once per segment into a
 *               CSR table (counting sort: offsets + entries, no linked lists).
 *               A query hashes its own k+1 segments and only looks at those
 *               buckets. Candidates are verified with one XOR and a popcount
 *               on a 2-bit-per-base code (20 letters = 40 bits).
 *               Duplicates (a site matching exactly in several segments) are
 *               removed without a hash set: a hit is reported only from the
 *               first segment that matches exactly.
 *
 * Build:  O(S * (k+1)) once, S = number of PAM sites (about n/8).
 * Query:  O(total bucket sizes), far below O(n).
 */
#include "crispr.h"
#include <stdlib.h>
#include <string.h>

#define M40 0x5555555555ULL /* low bit of each of the 20 two-bit slots */

int guide_encode(const char *s, uint64_t *code) {
    uint64_t c = 0;
    for (int j = 0; j < GUIDE_LEN; j++) {
        int b = base2(s[j]);
        if (b < 0) return -1;
        c |= (uint64_t)b << (2 * j);
    }
    if (s[GUIDE_LEN] != 0) return -1;
    *code = c;
    return 0;
}

void guide_decode(uint64_t code, char out[GUIDE_LEN + 1]) {
    static const char LET[4] = {'A', 'C', 'G', 'T'};
    for (int j = 0; j < GUIDE_LEN; j++) out[j] = LET[(code >> (2 * j)) & 3];
    out[GUIDE_LEN] = 0;
}

uint64_t site_code(const char *seq, size_t p, char strand) {
    uint64_t c = 0;
    if (strand == '+') {
        for (int j = 0; j < GUIDE_LEN; j++) c |= (uint64_t)base2(seq[p + (size_t)j]) << (2 * j);
    } else { /* reverse complement of seq[p+3 .. p+22]; complement = 3 - code */
        for (int j = 0; j < GUIDE_LEN; j++)
            c |= (uint64_t)(3 - base2(seq[p + 22 - (size_t)j])) << (2 * j);
    }
    return c;
}

static uint32_t mask20(uint64_t m) {
    uint32_t r = 0;
    for (int j = 0; j < GUIDE_LEN; j++) r |= (uint32_t)((m >> (2 * j)) & 1) << j;
    return r;
}

static uint64_t seg_mask(int len) { return (1ULL << (2 * len)) - 1; }
static uint64_t seg_key(uint64_t code, int start, int len) {
    return (code >> (2 * start)) & seg_mask(len);
}
static size_t bucket_of(uint64_t key, int seg, size_t nb) {
    return (size_t)(mix64(key * 8 + (uint64_t)seg) & (nb - 1));
}

int ot_build(const Genome *g, int k, OtIndex *ix) {
    memset(ix, 0, sizeof *ix);
    if (k < 0 || k > 5) return -1;
    ix->k = k;
    ix->nseg = k + 1;
    for (int i = 0; i < ix->nseg; i++) {
        ix->seg_start[i] = i * GUIDE_LEN / ix->nseg;
        ix->seg_len[i] = (i + 1) * GUIDE_LEN / ix->nseg - ix->seg_start[i];
    }
    size_t cap = 1024, ns = 0;
    Site *s = xmalloc(cap * sizeof *s);
    for (size_t p = 0; p + SITE_LEN <= g->n; p++) {
        if (g->seq[p + 21] == 'G' && g->seq[p + 22] == 'G') {
            if (ns == cap) { cap *= 2; s = realloc(s, cap * sizeof *s); if (!s) exit(1); }
            s[ns++] = (Site){p, '+', site_code(g->seq, p, '+')};
        }
        if (g->seq[p] == 'C' && g->seq[p + 1] == 'C') {
            if (ns == cap) { cap *= 2; s = realloc(s, cap * sizeof *s); if (!s) exit(1); }
            s[ns++] = (Site){p, '-', site_code(g->seq, p, '-')};
        }
    }
    ix->sites = s;
    ix->nsites = ns;
    size_t nb = 16;
    while (nb < ns) nb <<= 1;
    ix->nbuckets = nb;
    uint32_t *cursor = xmalloc(nb * sizeof *cursor);
    for (int i = 0; i < ix->nseg; i++) {
        uint32_t *off = xcalloc(nb + 1, sizeof *off);
        uint32_t *ent = xmalloc((ns ? ns : 1) * sizeof *ent);
        for (size_t t = 0; t < ns; t++)
            off[bucket_of(seg_key(s[t].code, ix->seg_start[i], ix->seg_len[i]), i, nb) + 1]++;
        for (size_t b = 0; b < nb; b++) off[b + 1] += off[b];
        memcpy(cursor, off, nb * sizeof *cursor);
        for (size_t t = 0; t < ns; t++) {
            size_t b = bucket_of(seg_key(s[t].code, ix->seg_start[i], ix->seg_len[i]), i, nb);
            ent[cursor[b]++] = (uint32_t)t;
        }
        ix->off[i] = off;
        ix->ent[i] = ent;
    }
    free(cursor);
    return 0;
}

void ot_free(OtIndex *ix) {
    for (int i = 0; i < 6; i++) { free(ix->off[i]); free(ix->ent[i]); }
    free(ix->sites);
    memset(ix, 0, sizeof *ix);
}

void ot_query(const OtIndex *ix, uint64_t guide, HitVec *out) {
    for (int i = 0; i < ix->nseg; i++) {
        uint64_t key = seg_key(guide, ix->seg_start[i], ix->seg_len[i]);
        size_t b = bucket_of(key, i, ix->nbuckets);
        for (uint32_t e = ix->off[i][b]; e < ix->off[i][b + 1]; e++) {
            const Site *st = &ix->sites[ix->ent[i][e]];
            uint64_t x = guide ^ st->code;
            uint64_t m = (x | (x >> 1)) & M40; /* one bit per mismatching letter */
            int mm = __builtin_popcountll(m);
            if (mm > ix->k) continue;
            /* segment i must match exactly (the bucket can hold hash collisions) */
            if ((m >> (2 * ix->seg_start[i])) & seg_mask(ix->seg_len[i])) continue;
            /* report only from the first exactly matching segment */
            int dup = 0;
            for (int j = 0; j < i && !dup; j++)
                if (((m >> (2 * ix->seg_start[j])) & seg_mask(ix->seg_len[j])) == 0) dup = 1;
            if (dup) continue;
            hitvec_push(out, (Hit){st->pos, st->strand, mm, mask20(m)});
        }
    }
}

void ot_brute(const Genome *g, uint64_t guide, int k, HitVec *out) {
    char gc[GUIDE_LEN + 1];
    guide_decode(guide, gc);
    for (size_t p = 0; p + SITE_LEN <= g->n; p++) {
        const char *w = g->seq + p;
        if (w[21] == 'G' && w[22] == 'G') {
            int mm = 0;
            uint32_t mask = 0;
            for (int j = 0; j < GUIDE_LEN && mm <= k; j++)
                if (w[j] != gc[j]) { mm++; mask |= 1u << j; }
            if (mm <= k) hitvec_push(out, (Hit){p, '+', mm, mask});
        }
        if (w[0] == 'C' && w[1] == 'C') {
            int mm = 0;
            uint32_t mask = 0;
            for (int j = 0; j < GUIDE_LEN && mm <= k; j++)
                if (complement(w[22 - j]) != gc[j]) { mm++; mask |= 1u << j; }
            if (mm <= k) hitvec_push(out, (Hit){p, '-', mm, mask});
        }
    }
}
