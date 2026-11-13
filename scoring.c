/* scoring.c: score every candidate guide in a region, keep the best K.
 *
 * Risk per off-target hit: published MIT (Hsu 2013) position weights, where a
 * mismatch near the PAM weighs more. For a hit with n >= 1 mismatches:
 *     risk = prod(1 - W[mismatch positions]) * 1/(4(19-d)/19 + 1) * 1/n^2
 * d = average gap between neighbouring mismatches (clustered mismatches are
 * more dangerous). A second perfect match elsewhere (n = 0) counts as risk 1.
 * specificity = 100 / (1 + sum of risks); GC outside 40-70% costs 20 points;
 * guides containing TTTT are rejected (a run of T stops the polymerase that
 * makes the guide).
 *
 * Heap top-K:          O(N log K), keeps a size-K heap whose root is the worst
 *                      kept guide.
 * Quicksort (final):   O(N log N) average, median-of-three pivot, insertion
 *                      sort on short ranges, recurse into the smaller half.
 * Median of medians:   O(N) worst case. Groups of 5, median of the group
 *                      medians as pivot, three-way partition.
 */
#include "crispr.h"
#include <stdlib.h>
#include <string.h>

typedef int (*RankCmp)(const Ranked *, const Ranked *);

static const double W[GUIDE_LEN] = {0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
                                    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,
                                    0.804, 0.685, 0.583};

double hit_risk(const Hit *h) {
    if (h->mm == 0) return 1.0;
    double prod = 1.0;
    int first = -1, last = -1;
    for (int j = 0; j < GUIDE_LEN; j++) {
        if (h->mask & (1u << j)) {
            prod *= 1.0 - W[j];
            if (first < 0) first = j;
            last = j;
        }
    }
    double d = h->mm >= 2 ? (double)(last - first) / (h->mm - 1) : 19.0;
    double dterm = 1.0 / ((19.0 - d) / 19.0 * 4.0 + 1.0);
    return prod * dterm / ((double)h->mm * h->mm);
}

/* best first: higher score, then lower position, then strand */
static int cmp_best(const Ranked *a, const Ranked *b) {
    if (a->score != b->score) return a->score > b->score ? -1 : 1;
    if (a->pos != b->pos) return a->pos < b->pos ? -1 : 1;
    if (a->strand != b->strand) return a->strand < b->strand ? -1 : 1;
    return 0;
}
static int cmp_worst(const Ranked *a, const Ranked *b) { return cmp_best(b, a); }

/* ---------- bounded heap: keeps the k smallest items by cmp ---------- */
static void swap_idx(size_t *a, size_t *b) { size_t t = *a; *a = *b; *b = t; }

static void sift_down(size_t *h, size_t n, size_t i, const Ranked *it, RankCmp cmp) {
    for (;;) {
        size_t l = 2 * i + 1, r = l + 1, big = i;
        if (l < n && cmp(&it[h[l]], &it[h[big]]) > 0) big = l;
        if (r < n && cmp(&it[h[r]], &it[h[big]]) > 0) big = r;
        if (big == i) return;
        swap_idx(&h[i], &h[big]);
        i = big;
    }
}

/* h must have room for k indices. Returns how many were kept. Root h[0] is the
 * largest (worst) kept item, so a new item only has to beat the root. */
static size_t topk_select(const Ranked *it, size_t n, size_t k, RankCmp cmp, size_t *h) {
    size_t sz = 0;
    if (k == 0) return 0;
    for (size_t i = 0; i < n; i++) {
        if (sz < k) {
            size_t c = sz++;
            h[c] = i;
            while (c > 0) {
                size_t p = (c - 1) / 2;
                if (cmp(&it[h[c]], &it[h[p]]) > 0) { swap_idx(&h[c], &h[p]); c = p; }
                else break;
            }
        } else if (cmp(&it[i], &it[h[0]]) < 0) {
            h[0] = i;
            sift_down(h, sz, 0, it, cmp);
        }
    }
    return sz;
}

/* ---------- quicksort on Ranked ---------- */
static void rswap(Ranked *a, Ranked *b) { Ranked t = *a; *a = *b; *b = t; }

static void quicksort_ranked(Ranked *a, long lo, long hi, RankCmp cmp) {
    while (hi - lo > 16) {
        long mid = lo + (hi - lo) / 2;
        if (cmp(&a[mid], &a[lo]) < 0) rswap(&a[mid], &a[lo]);
        if (cmp(&a[hi], &a[lo]) < 0) rswap(&a[hi], &a[lo]);
        if (cmp(&a[hi], &a[mid]) < 0) rswap(&a[hi], &a[mid]);
        Ranked pivot = a[mid];
        long i = lo, j = hi;
        while (i <= j) {
            while (cmp(&a[i], &pivot) < 0) i++;
            while (cmp(&a[j], &pivot) > 0) j--;
            if (i <= j) { rswap(&a[i], &a[j]); i++; j--; }
        }
        if (j - lo < hi - i) { quicksort_ranked(a, lo, j, cmp); lo = i; }
        else { quicksort_ranked(a, i, hi, cmp); hi = j; }
    }
    for (long i = lo + 1; i <= hi; i++) { /* insertion sort */
        Ranked key = a[i];
        long j = i - 1;
        while (j >= lo && cmp(&key, &a[j]) < 0) { a[j + 1] = a[j]; j--; }
        a[j + 1] = key;
    }
}

/* ---------- median of medians (BFPRT) ---------- */
static void dswap(double *a, double *b) { double t = *a; *a = *b; *b = t; }
static void dinsertion(double *a, size_t lo, size_t hi) { /* sorts [lo, hi) */
    for (size_t i = lo + 1; i < hi; i++) {
        double key = a[i];
        size_t j = i;
        while (j > lo && a[j - 1] > key) { a[j] = a[j - 1]; j--; }
        a[j] = key;
    }
}

/* Returns the element that would sit at index t if a[0..n) were sorted.
 * Rearranges a. Worst case O(n). */
static double mom_select(double *a, size_t n, size_t t) {
    size_t lo = 0, hi = n;
    for (;;) {
        size_t m = hi - lo;
        if (m <= 5) { dinsertion(a, lo, hi); return a[t]; }
        size_t ng = (m + 4) / 5;
        for (size_t g = 0; g < ng; g++) {
            size_t s = lo + g * 5, e = s + 5 < hi ? s + 5 : hi;
            dinsertion(a, s, e);
            dswap(&a[lo + g], &a[s + (e - s - 1) / 2]); /* group median to the front */
        }
        double pivot = mom_select(a + lo, ng, (ng - 1) / 2);
        size_t lt = lo, i = lo, gt = hi; /* <pivot | ==pivot | >pivot */
        while (i < gt) {
            if (a[i] < pivot) dswap(&a[lt++], &a[i++]);
            else if (a[i] > pivot) dswap(&a[i], &a[--gt]);
            else i++;
        }
        if (t < lt) hi = lt;
        else if (t < gt) return pivot;
        else lo = gt;
    }
}

static double median_of(const double *v, size_t n) {
    if (n == 0) return 0.0;
    double *c = xmalloc(n * sizeof *c);
    memcpy(c, v, n * sizeof *c);
    double med = mom_select(c, n, n / 2);
    if (n % 2 == 0) {
        memcpy(c, v, n * sizeof *c);
        med = (med + mom_select(c, n, n / 2 - 1)) / 2.0;
    }
    free(c);
    return med;
}

/* ---------- ranking ---------- */
static int has_tttt(const char *g) { return strstr(g, "TTTT") != NULL; }

int rank_region(const Genome *g, const OtIndex *ix, size_t from, size_t to,
                size_t K, RankResult *res) {
    memset(res, 0, sizeof *res);
    if (to > g->n) to = g->n;
    if (from >= to || to - from < SITE_LEN || K == 0) return -1;

    size_t cap = 256, ns = 0;
    Ranked *items = xmalloc(cap * sizeof *items);
    HitVec hv = {0};

    for (size_t p = from; p + SITE_LEN <= to; p++) {
        const char *w = g->seq + p;
        for (int s = 0; s < 2; s++) {
            char strand = s == 0 ? '+' : '-';
            if (strand == '+' ? !(w[21] == 'G' && w[22] == 'G') : !(w[0] == 'C' && w[1] == 'C'))
                continue;
            res->ncand++;
            Ranked r;
            memset(&r, 0, sizeof r);
            r.pos = p;
            r.strand = strand;
            if (strand == '+') {
                memcpy(r.guide, w, GUIDE_LEN);
                memcpy(r.pam, w + 20, 3);
            } else {
                for (int j = 0; j < GUIDE_LEN; j++) r.guide[j] = complement(w[22 - j]);
                r.pam[0] = complement(w[2]); r.pam[1] = complement(w[1]); r.pam[2] = complement(w[0]);
            }
            r.guide[GUIDE_LEN] = 0;
            r.pam[3] = 0;
            if (has_tttt(r.guide)) { res->nrej_tttt++; continue; }
            int gcn = 0;
            for (int j = 0; j < GUIDE_LEN; j++) gcn += (r.guide[j] == 'G' || r.guide[j] == 'C');
            r.gc = gcn / (double)GUIDE_LEN;

            uint64_t code = site_code(g->seq, p, strand);
            hv.n = 0;
            ot_query(ix, code, &hv);
            for (size_t h = 0; h < hv.n; h++) {
                if (hv.a[h].pos == p && hv.a[h].strand == strand) continue; /* the target itself */
                r.risk += hit_risk(&hv.a[h]);
                r.nhits++;
            }
            r.spec = 100.0 / (1.0 + r.risk);
            r.score = r.spec - ((r.gc < 0.40 || r.gc > 0.70) ? 20.0 : 0.0);
            if (r.score < 0) r.score = 0;
            if (ns == cap) { cap *= 2; items = realloc(items, cap * sizeof *items); if (!items) exit(1); }
            items[ns++] = r;
        }
    }
    hitvec_free(&hv);
    res->nscored = ns;
    if (ns == 0) { free(items); return 0; }

    size_t *heap = xmalloc((K < ns ? K : ns) * sizeof *heap);
    size_t kept = topk_select(items, ns, K, cmp_best, heap);
    res->top = xmalloc(kept * sizeof *res->top);
    for (size_t i = 0; i < kept; i++) res->top[i] = items[heap[i]];
    res->ntop = kept;
    quicksort_ranked(res->top, 0, (long)kept - 1, cmp_best);
    free(heap);

    size_t h3[3];
    res->nworst = topk_select(items, ns, 3, cmp_worst, h3);
    for (size_t i = 0; i < res->nworst; i++) res->worst[i] = items[h3[i]];
    quicksort_ranked(res->worst, 0, (long)res->nworst - 1, cmp_worst);

    double *specs = xmalloc(ns * sizeof *specs);
    for (size_t i = 0; i < ns; i++) specs[i] = items[i].spec;
    res->median_spec = median_of(specs, ns);
    free(specs);
    free(items);
    return 0;
}

void rank_free(RankResult *res) { free(res->top); memset(res, 0, sizeof *res); }
