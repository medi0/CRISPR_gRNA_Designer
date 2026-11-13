/* Temporary Phase 1 test driver. Replaced by `crispr demo` in Phase 2. */
#include "../crispr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LLU(x) ((unsigned long long)(x))

static int fails = 0;
#define CHECK(c, msg) do { if (!(c)) { printf("  FAIL: %s\n", msg); fails++; } else printf("  ok:   %s\n", msg); } while (0)

static void test_matchers(const Genome *g, const SynthInfo *si) {
    MatchFn fn[3] = {kmp_search, bm_search, rk_search};
    const char *nm[3] = {"KMP", "Boyer-Moore", "Rabin-Karp"};
    size_t cnt[3];
    uint64_t sig[3];
    printf("[target finder] genome %llu bases\n", LLU(g->n));
    for (int i = 0; i < 3; i++) {
        TargetVec tv = {0};
        double t0 = now_ms();
        find_targets(g, fn[i], &tv);
        double ms = now_ms() - t0;
        cnt[i] = tv.n;
        sig[i] = targetvec_signature(&tv);
        printf("  %-12s targets=%llu sig=%016llx  %.2f ms\n", nm[i], LLU(tv.n), LLU(sig[i]), ms);
        targetvec_free(&tv);
    }
    CHECK(cnt[0] == cnt[1] && cnt[1] == cnt[2] && sig[0] == sig[1] && sig[1] == sig[2],
          "all three matchers: same count and signature");

    char pat[17];
    memcpy(pat, si->guide, 16);
    pat[16] = 0;
    size_t pc[3];
    printf("  16-letter pattern %s:\n", pat);
    for (int i = 0; i < 3; i++) {
        double t0 = now_ms();
        size_t c = 0;
        for (int rep = 0; rep < 5; rep++) c = fn[i](g->seq, g->n, pat, 16, NULL);
        printf("    %-12s matches=%llu  %.2f ms per scan\n", nm[i], LLU(c), (now_ms() - t0) / 5);
        pc[i] = c;
    }
    CHECK(pc[0] == pc[1] && pc[1] == pc[2] && pc[0] >= 1, "16-letter pattern: matchers agree");
}

static int same_hits(HitVec *a, HitVec *b) {
    hitvec_sort(a);
    hitvec_sort(b);
    if (a->n != b->n || hitvec_signature(a) != hitvec_signature(b)) return 0;
    for (size_t i = 0; i < a->n; i++)
        if (a->a[i].pos != b->a[i].pos || a->a[i].strand != b->a[i].strand ||
            a->a[i].mm != b->a[i].mm || a->a[i].mask != b->a[i].mask) return 0;
    return 1;
}

static void test_offtarget(const Genome *g, const SynthInfo *si) {
    printf("[off-target] genome %llu bases\n", LLU(g->n));
    OtIndex ix;
    double t0 = now_ms();
    ot_build(g, 3, &ix);
    double build = now_ms() - t0;
    printf("  index k=3: %llu sites, %d segments, build %.1f ms\n", LLU(ix.nsites), ix.nseg, build);

    TargetVec tv = {0};
    find_targets(g, kmp_search, &tv);
    Rng r;
    rng_seed(&r, 777);
    const int Q = 20;
    uint64_t codes[20];
    guide_encode(si->guide, &codes[0]);
    for (int q = 1; q < Q; q++) {
        Target t = tv.a[rng_below(&r, (uint32_t)tv.n)];
        codes[q] = site_code(g->seq, t.pos, t.strand);
    }
    double tb = 0, ti = 0;
    int equal = 0;
    for (int q = 0; q < Q; q++) {
        HitVec a = {0}, b = {0};
        double s0 = now_ms();
        ot_brute(g, codes[q], 3, &a);
        double s1 = now_ms();
        ot_query(&ix, codes[q], &b);
        double s2 = now_ms();
        tb += s1 - s0;
        ti += s2 - s1;
        if (same_hits(&a, &b)) equal++;
        else printf("  mismatch on guide %d: brute %llu index %llu\n", q, LLU(a.n), LLU(b.n));
        if (q == 0) {
            int by_mm[4] = {0};
            for (size_t i = 0; i < b.n; i++) by_mm[b.a[i].mm]++;
            printf("  planted guide at k=3: %llu hits (mm0=%d mm1=%d mm2=%d mm3=%d)\n",
                   LLU(b.n), by_mm[0], by_mm[1], by_mm[2], by_mm[3]);
            CHECK(b.n == 13 && by_mm[0] == 1, "planted guide: 1 intended site + 12 decoys at k=3");
            int found = 0;
            for (int d = 0; d < 12; d++)
                for (size_t i = 0; i < b.n; i++)
                    if (b.a[i].pos == si->decoy_pos[d] && b.a[i].strand == si->decoy_strand[d] &&
                        b.a[i].mm == si->decoy_mm[d]) { found++; break; }
            CHECK(found == 12, "every planted decoy found at its planted position and mismatch count");
        }
        hitvec_free(&a);
        hitvec_free(&b);
    }
    printf("  identical results: %d of %d guides\n", equal, Q);
    CHECK(equal == Q, "20 guides: index equals brute force (count, signature, every hit)");
    printf("  brute %.2f ms/query, index %.4f ms/query, speed-up %.0fx per query\n", tb / Q, ti / Q, tb / ti);
    printf("  including build, %d queries: brute %.1f ms vs index %.1f ms -> %.1fx; break-even at %.1f queries\n",
           Q, tb, build + ti, tb / (build + ti), build / ((tb - ti) / Q));
    ot_free(&ix);

    /* other k values on the planted guide and a few random guides */
    int allk = 1;
    for (int k = 0; k <= 5; k++) {
        OtIndex ik;
        ot_build(g, k, &ik);
        for (int q = 0; q < 4; q++) {
            HitVec a = {0}, b = {0};
            ot_brute(g, codes[q], k, &a);
            ot_query(&ik, codes[q], &b);
            if (!same_hits(&a, &b)) allk = 0;
            hitvec_free(&a);
            hitvec_free(&b);
        }
        ot_free(&ik);
    }
    CHECK(allk, "k = 0..5, 4 guides each: index equals brute force");
    targetvec_free(&tv);
}

static void test_rank(const Genome *g, const SynthInfo *si) {
    printf("[ranking]\n");
    OtIndex ix;
    ot_build(g, 3, &ix);
    size_t from = si->guide_pos - 15000, to = si->guide_pos + 15000;
    RankResult rr;
    double t0 = now_ms();
    int rc = rank_region(g, &ix, from, to, 10, &rr);
    printf("  region %llu..%llu: %llu candidates, %llu scored, %llu rejected (TTTT), %.1f ms, median specificity %.2f\n",
           LLU(from), LLU(to), LLU(rr.ncand), LLU(rr.nscored), LLU(rr.nrej_tttt), now_ms() - t0, rr.median_spec);
    CHECK(rc == 0 && rr.ntop == 10, "top 10 returned");
    for (size_t i = 0; i < 5 && i < rr.ntop; i++)
        printf("  top%llu  %c %llu %s %s score=%.2f hits=%d\n", LLU(i + 1), rr.top[i].strand,
               LLU(rr.top[i].pos), rr.top[i].guide, rr.top[i].pam, rr.top[i].score, rr.top[i].nhits);
    int planted_worst = 0;
    for (size_t i = 0; i < rr.nworst; i++) {
        printf("  worst%llu %c %llu %s %s spec=%.2f score=%.2f hits=%d\n", LLU(i + 1), rr.worst[i].strand,
               LLU(rr.worst[i].pos), rr.worst[i].guide, rr.worst[i].pam, rr.worst[i].spec,
               rr.worst[i].score, rr.worst[i].nhits);
        if (rr.worst[i].pos == si->guide_pos && rr.worst[i].strand == '+') planted_worst = 1;
    }
    CHECK(planted_worst, "planted guide is among the 3 worst");
    CHECK(rr.nworst == 3 && rr.worst[0].score <= rr.worst[1].score && rr.worst[1].score <= rr.worst[2].score,
          "worst list sorted worst first");
    int sorted = 1, top100 = 0;
    for (size_t i = 1; i < rr.ntop; i++) if (rr.top[i - 1].score < rr.top[i].score) sorted = 0;
    for (size_t i = 0; i < rr.ntop; i++) if (rr.top[i].score == 100.0 && rr.top[i].nhits == 0) top100++;
    CHECK(sorted, "top list sorted best first (quicksort)");
    CHECK(top100 == 10, "guides with no near-matches score exactly 100");

    /* cross-check heap + quicksort against a plain full sort of all scores */
    RankResult all;
    rank_region(g, &ix, from, to, rr.nscored, &all);
    int ok = all.ntop == rr.nscored;
    for (size_t i = 0; ok && i < rr.ntop; i++) if (all.top[i].pos != rr.top[i].pos) ok = 0;
    CHECK(ok, "heap top-10 equals first 10 of the fully sorted list");
    double *s = malloc(all.ntop * sizeof *s);
    for (size_t i = 0; i < all.ntop; i++) s[i] = all.top[i].spec;
    /* median from the sorted-by-score list is not sorted by spec; sort a copy */
    for (size_t i = 1; i < all.ntop; i++) { double k = s[i]; size_t j = i; while (j && s[j - 1] > k) { s[j] = s[j - 1]; j--; } s[j] = k; }
    size_t n = all.ntop;
    double med = n % 2 ? s[n / 2] : (s[n / 2 - 1] + s[n / 2]) / 2;
    printf("  median check: median-of-medians %.6f, plain sort %.6f\n", rr.median_spec, med);
    CHECK(rr.median_spec == med, "median of medians equals median from a full sort");
    free(s);
    rank_free(&all);
    rank_free(&rr);
    ot_free(&ix);
}

int main(void) {
    Genome g;
    SynthInfo si;
    /* FASTA round trip and error path */
    char err[200];
    const char *messy = ">test seq\r\nacgtNNACGT\r\nTTGA\r\n>second\nGGGG\n";
    CHECK(genome_parse(messy, strlen(messy), &g, err, sizeof err) == 0 && g.n == 12 &&
          strcmp(g.seq, "ACGTACGTTTGA") == 0, "FASTA: first record only, upper-cased, N and CR dropped");
    genome_free(&g);
    CHECK(genome_parse("hello world", 11, &g, err, sizeof err) != 0, "garbage input gives an error");
    printf("  error text: %s\n", err);

    CHECK(synth_genome(300000, 42, &g, &si) == 0, "synthetic 300 kb genome");
    test_matchers(&g, &si);
    test_offtarget(&g, &si);
    genome_free(&g);

    CHECK(synth_genome(4600000, 42, &g, &si) == 0, "synthetic 4.6 Mb genome");
    test_matchers(&g, &si);
    test_offtarget(&g, &si);
    test_rank(&g, &si);
    genome_free(&g);
    printf("\n%s (%d failures)\n", fails ? "SOME TESTS FAILED" : "ALL PHASE 1 TESTS PASSED", fails);
    return fails != 0;
}
