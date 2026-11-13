/* align.c: sequence alignment with match +2, mismatch -1, gap -2.
 *
 * Needleman-Wunsch (global): fill a (m+1) x (n+1) table, trace back from the
 *     corner. O(m*n) time and O(m*n) space (16-bit cells here, so the table
 *     is refused when m+n > 16000, where scores could overflow).
 * Hirschberg: same optimal score, linear space. Run the NW recurrence keeping
 *     one row, forwards over the first half of A and backwards over the second
 *     half. The best split point of B is where forward + backward scores add up
 *     highest. Recurse on both halves. O(m*n) time (about twice NW), O(n) space.
 * Smith-Waterman (local): scores can not go below 0; answer is the best cell.
 *     Score only, O(m*n) time, O(n) space.
 */
#include "crispr.h"
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#endif

size_t peak_rss_bytes(void) {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc)) return (size_t)pmc.PeakWorkingSetSize;
    return 0;
#else
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (size_t)ru.ru_maxrss * 1024; /* Linux reports kilobytes */
#endif
}

static int sc(char x, char y) { return x == y ? AL_MATCH : AL_MISMATCH; }
static int max3(int a, int b, int c) { int m = a > b ? a : b; return m > c ? m : c; }

void alignment_free(Alignment *al) { free(al->a); free(al->b); memset(al, 0, sizeof *al); }

int alignment_score(const Alignment *al) {
    int s = 0;
    for (size_t i = 0; i < al->len; i++) {
        if (al->a[i] == '-' || al->b[i] == '-') s += AL_GAP;
        else s += sc(al->a[i], al->b[i]);
    }
    return s;
}

static void reverse_str(char *s, size_t n) {
    for (size_t i = 0; i < n / 2; i++) { char t = s[i]; s[i] = s[n - 1 - i]; s[n - 1 - i] = t; }
}

int nw_align(const char *a, size_t m, const char *b, size_t n, Alignment *out) {
    memset(out, 0, sizeof *out);
    if (m + n > 16000) return -1;
    size_t W = n + 1;
    int16_t *T = xmalloc((m + 1) * W * sizeof *T);
    for (size_t j = 0; j <= n; j++) T[j] = (int16_t)(AL_GAP * (int)j);
    for (size_t i = 0; i <= m; i++) T[i * W] = (int16_t)(AL_GAP * (int)i);
    for (size_t i = 1; i <= m; i++)
        for (size_t j = 1; j <= n; j++)
            T[i * W + j] = (int16_t)max3(T[(i - 1) * W + j - 1] + sc(a[i - 1], b[j - 1]),
                                        T[(i - 1) * W + j] + AL_GAP, T[i * W + j - 1] + AL_GAP);
    char *ra = xmalloc(m + n + 1), *rb = xmalloc(m + n + 1);
    size_t i = m, j = n, k = 0;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && T[i * W + j] == T[(i - 1) * W + j - 1] + sc(a[i - 1], b[j - 1])) {
            ra[k] = a[i - 1]; rb[k] = b[j - 1]; i--; j--;
        } else if (i > 0 && T[i * W + j] == T[(i - 1) * W + j] + AL_GAP) {
            ra[k] = a[i - 1]; rb[k] = '-'; i--;
        } else {
            ra[k] = '-'; rb[k] = b[j - 1]; j--;
        }
        k++;
    }
    ra[k] = rb[k] = 0;
    reverse_str(ra, k);
    reverse_str(rb, k);
    out->score = T[m * W + n];
    out->a = ra; out->b = rb; out->len = k;
    out->bytes = (m + 1) * W * sizeof *T + 2 * (m + n + 1);
    free(T);
    return 0;
}

/* one-row NW: row[j] = best score of aligning A (or reversed A) against the
 * first j letters of B (or reversed B) */
static void last_row(const char *a, size_t m, const char *b, size_t n, int rev, int *row) {
    for (size_t j = 0; j <= n; j++) row[j] = AL_GAP * (int)j;
    for (size_t i = 1; i <= m; i++) {
        int diag = row[0];
        row[0] = AL_GAP * (int)i;
        char ai = rev ? a[m - i] : a[i - 1];
        for (size_t j = 1; j <= n; j++) {
            char bj = rev ? b[n - j] : b[j - 1];
            int up = row[j];
            row[j] = max3(diag + sc(ai, bj), up + AL_GAP, row[j - 1] + AL_GAP);
            diag = up;
        }
    }
}

typedef struct { int *L, *R; char *oa, *ob; size_t k; } Hb;

static void hb_put(Hb *h, char x, char y) { h->oa[h->k] = x; h->ob[h->k] = y; h->k++; }

static void hb_rec(Hb *h, const char *a, size_t m, const char *b, size_t n) {
    if (m == 0) { for (size_t j = 0; j < n; j++) hb_put(h, '-', b[j]); return; }
    if (n == 0) { for (size_t i = 0; i < m; i++) hb_put(h, a[i], '-'); return; }
    if (m == 1) {
        /* one letter of A: pairing it with some letter of B always beats two gaps,
         * so pick the first matching letter (or the first letter if none matches) */
        size_t at = 0;
        for (size_t j = 0; j < n; j++) if (b[j] == a[0]) { at = j; break; }
        for (size_t j = 0; j < n; j++) {
            if (j == at) hb_put(h, a[0], b[j]);
            else hb_put(h, '-', b[j]);
        }
        return;
    }
    size_t mid = m / 2;
    last_row(a, mid, b, n, 0, h->L);
    last_row(a + mid, m - mid, b, n, 1, h->R);
    size_t best = 0;
    int bs = h->L[0] + h->R[n];
    for (size_t j = 1; j <= n; j++) {
        int v = h->L[j] + h->R[n - j];
        if (v > bs) { bs = v; best = j; }
    }
    hb_rec(h, a, mid, b, best);
    hb_rec(h, a + mid, m - mid, b + best, n - best);
}

int hirschberg_align(const char *a, size_t m, const char *b, size_t n, Alignment *out) {
    memset(out, 0, sizeof *out);
    Hb h;
    h.L = xmalloc((n + 1) * sizeof *h.L);
    h.R = xmalloc((n + 1) * sizeof *h.R);
    h.oa = xmalloc(m + n + 1);
    h.ob = xmalloc(m + n + 1);
    h.k = 0;
    hb_rec(&h, a, m, b, n);
    h.oa[h.k] = h.ob[h.k] = 0;
    out->a = h.oa; out->b = h.ob; out->len = h.k;
    out->score = alignment_score(out);
    out->bytes = 2 * (n + 1) * sizeof *h.L + 2 * (m + n + 1);
    free(h.L); free(h.R);
    return 0;
}

int sw_score(const char *a, size_t m, const char *b, size_t n, size_t *bytes) {
    int *row = xcalloc(n + 1, sizeof *row);
    int best = 0;
    for (size_t i = 1; i <= m; i++) {
        int diag = 0; /* row[0] stays 0 */
        for (size_t j = 1; j <= n; j++) {
            int up = row[j];
            int v = max3(diag + sc(a[i - 1], b[j - 1]), up + AL_GAP, row[j - 1] + AL_GAP);
            if (v < 0) v = 0;
            row[j] = v;
            if (v > best) best = v;
            diag = up;
        }
    }
    if (bytes) *bytes = (n + 1) * sizeof *row;
    free(row);
    return best;
}
