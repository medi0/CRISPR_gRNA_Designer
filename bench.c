/* bench.c: timing runs behind the `bench` command and /api/bench.
 * Measures wall-clock time with now_ms(); returns rows, prints nothing. */
#include "crispr.h"
#include <stdlib.h>

int bench_genomes(const size_t *sizes, int ns, GenomeBenchRow *rows) {
    for (int i = 0; i < ns; i++) {
        Genome g;
        SynthInfo si;
        GenomeBenchRow *r = &rows[i];
        if (synth_genome(sizes[i], 42, &g, &si) != 0) return -1;
        r->n = g.n;
        MatchFn fn[3] = {kmp_search, bm_search, rk_search};
        double *slot[3] = {&r->kmp_ms, &r->bm_ms, &r->rk_ms};
        TargetVec keep = {0};
        for (int m = 0; m < 3; m++) {
            TargetVec tv = {0};
            double t0 = now_ms();
            find_targets(&g, fn[m], &tv);
            *slot[m] = now_ms() - t0;
            if (m == 0) keep = tv; else targetvec_free(&tv);
        }
        r->sites = keep.n;
        IslandVec iv = {0};
        double t0 = now_ms();
        cpg_islands(g.seq, g.n, &iv);
        r->cpg_ms = now_ms() - t0;
        islandvec_free(&iv);
        OtIndex ix;
        t0 = now_ms();
        ot_build(&g, 3, &ix);
        r->build_ms = now_ms() - t0;
        Rng rg;
        rng_seed(&rg, 99);
        const int Q = 5;
        double tb = 0, tq = 0;
        for (int q = 0; q < Q; q++) {
            Target t = keep.a[rng_below(&rg, (uint32_t)keep.n)];
            uint64_t code = site_code(g.seq, t.pos, t.strand);
            HitVec a = {0}, b = {0};
            double s0 = now_ms();
            ot_brute(&g, code, 3, &a);
            double s1 = now_ms();
            ot_query(&ix, code, &b);
            double s2 = now_ms();
            tb += s1 - s0; tq += s2 - s1;
            hitvec_free(&a); hitvec_free(&b);
        }
        r->brute_ms = tb / Q;
        r->query_ms = tq / Q;
        ot_free(&ix);
        targetvec_free(&keep);
        genome_free(&g);
    }
    return 0;
}

int bench_align(const size_t *lens, int nl, AlignBenchRow *rows) {
    static const char L[4] = {'A', 'C', 'G', 'T'};
    Rng r;
    rng_seed(&r, 5);
    for (int i = 0; i < nl; i++) {
        size_t n = lens[i];
        if (n > 8000) return -1;
        char *a = xmalloc(n + 1), *b = xmalloc(n + 1);
        for (size_t j = 0; j < n; j++) { a[j] = L[rng_below(&r, 4)]; b[j] = L[rng_below(&r, 4)]; }
        a[n] = b[n] = 0;
        Alignment x, y;
        double t0 = now_ms();
        hirschberg_align(a, n, b, n, &y);
        double t1 = now_ms();
        nw_align(a, n, b, n, &x);
        double t2 = now_ms();
        rows[i] = (AlignBenchRow){n, t2 - t1, t1 - t0, x.bytes, y.bytes};
        alignment_free(&x); alignment_free(&y);
        free(a); free(b);
    }
    return 0;
}
