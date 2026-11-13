/* strmatch.c: three exact string matchers and the CRISPR target finder.
 *
 * KMP          O(n + m) always. Builds a failure table so the text pointer
 *              never moves backwards.
 * Boyer-Moore  bad-character rule only. Best case O(n/m), worst case O(n*m).
 *              Compares right to left and jumps over text it can rule out.
 * Rabin-Karp   rolling hash. O(n + m) expected, O(n*m) if hashes collide a lot.
 *              Every hash hit is verified with memcmp.
 *
 * Target finder: a PAM-adjacent site is found by searching for the fixed part
 * of the PAM. "GG" on the forward strand, "CC" on the reverse strand (the
 * reverse strand reads CCN). The N and the 20 letters come for free from the
 * position. Cost: one pass of the chosen matcher per strand.
 */
#include "crispr.h"
#include <stdlib.h>
#include <string.h>

size_t kmp_search(const char *t, size_t n, const char *p, size_t m, PosVec *out) {
    if (m == 0 || m > n) return 0;
    size_t *fail = xmalloc(m * sizeof *fail);
    fail[0] = 0;
    for (size_t i = 1, k = 0; i < m; i++) {
        while (k > 0 && p[i] != p[k]) k = fail[k - 1];
        if (p[i] == p[k]) k++;
        fail[i] = k;
    }
    size_t cnt = 0;
    for (size_t i = 0, k = 0; i < n; i++) {
        while (k > 0 && t[i] != p[k]) k = fail[k - 1];
        if (t[i] == p[k]) k++;
        if (k == m) {
            cnt++;
            if (out) posvec_push(out, i + 1 - m);
            k = fail[k - 1];
        }
    }
    free(fail);
    return cnt;
}

size_t bm_search(const char *t, size_t n, const char *p, size_t m, PosVec *out) {
    if (m == 0 || m > n) return 0;
    long last[256];
    for (int c = 0; c < 256; c++) last[c] = -1;
    for (size_t j = 0; j < m; j++) last[(unsigned char)p[j]] = (long)j;
    size_t cnt = 0, s = 0;
    while (s <= n - m) {
        long j = (long)m - 1;
        while (j >= 0 && p[j] == t[s + (size_t)j]) j--;
        if (j < 0) {
            cnt++;
            if (out) posvec_push(out, s);
            s++;
        } else {
            long shift = j - last[(unsigned char)t[s + (size_t)j]];
            s += shift > 1 ? (size_t)shift : 1;
        }
    }
    return cnt;
}

size_t rk_search(const char *t, size_t n, const char *p, size_t m, PosVec *out) {
    if (m == 0 || m > n) return 0;
    const uint64_t B = 131, M = 1000000007ULL;
    uint64_t hp = 0, ht = 0, pw = 1; /* pw = B^(m-1) mod M */
    for (size_t i = 0; i < m; i++) {
        hp = (hp * B + (unsigned char)p[i]) % M;
        ht = (ht * B + (unsigned char)t[i]) % M;
        if (i + 1 < m) pw = pw * B % M;
    }
    size_t cnt = 0;
    for (size_t s = 0;; s++) {
        if (hp == ht && memcmp(t + s, p, m) == 0) {
            cnt++;
            if (out) posvec_push(out, s);
        }
        if (s + m >= n) break;
        ht = (ht + M - (uint64_t)(unsigned char)t[s] * pw % M) % M;
        ht = (ht * B + (unsigned char)t[s + m]) % M;
    }
    return cnt;
}

void find_targets(const Genome *g, MatchFn fn, TargetVec *out) {
    PosVec gg = {0}, cc = {0};
    fn(g->seq, g->n, "GG", 2, &gg);
    fn(g->seq, g->n, "CC", 2, &cc);
    for (size_t i = 0; i < gg.n; i++) {
        size_t at = gg.a[i];              /* "GG" starts at window pos + 21 */
        if (at >= 21) targetvec_push(out, (Target){at - 21, '+'});
    }
    for (size_t i = 0; i < cc.n; i++) {
        size_t p = cc.a[i];               /* "CC" starts the window */
        if (p + SITE_LEN <= g->n) targetvec_push(out, (Target){p, '-'});
    }
    posvec_free(&gg);
    posvec_free(&cc);
}
