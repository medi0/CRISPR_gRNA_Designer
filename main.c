/* main.c: command-line interface. All printing for the CLI happens here. */
#include "crispr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LLU(x) ((unsigned long long)(x))

static int checks = 0, check_fail = 0;
static void check(int ok, const char *msg) {
    checks++;
    if (!ok) check_fail++;
    printf("  [%s] %s\n", ok ? "PASS" : "FAIL", msg);
}

static int fail(const char *msg) { fprintf(stderr, "error: %s\n", msg); return 1; }

static int parse_ull(const char *s, unsigned long long *v) {
    char *end;
    if (!s || !*s) return -1;
    *v = strtoull(s, &end, 10);
    return *end ? -1 : 0;
}

/* "synth:N" makes a synthetic genome, anything else is a FASTA path */
static int load_genome_arg(const char *arg, Genome *g) {
    char err[256];
    if (!strncmp(arg, "synth:", 6)) {
        unsigned long long n;
        if (parse_ull(arg + 6, &n) != 0 || synth_genome((size_t)n, 42, g, NULL) != 0) {
            fprintf(stderr, "error: synth:N needs N of at least 20000.\n");
            return -1;
        }
        return 0;
    }
    if (genome_load(arg, g, err, sizeof err) != 0) { fprintf(stderr, "error: %s\n", err); return -1; }
    return 0;
}

static void usage(void) {
    puts("CRISPR Guide RNA Designer with a Secure Genome Vault\n"
         "\n"
         "usage: crispr <command> [arguments]\n"
         "\n"
         "A <genome> is a FASTA file, or synth:N for a built-in synthetic genome of N bases.\n"
         "\n"
         "  demo [n]                          run every module on a synthetic genome (default 300000)\n"
         "  gen <n> <out.fa> [seed]           write a synthetic genome as FASTA\n"
         "  find <genome>                     find NGG targets with KMP, Boyer-Moore, Rabin-Karp\n"
         "  offtarget <genome> <guide> [k]    every site within k mismatches (0 to 5, default 3)\n"
         "  rank <genome> <from> <to> [K] [k] best guides in a region (default K=10, k=3)\n"
         "  cpg <genome>                      CpG islands\n"
         "  lcs <genomeA> <genomeB> [len]     longest common subsequence of the first len letters\n"
         "  align <A> <B>                     A and B are letters or rand:LENGTH\n"
         "  keygen <bits> <prefix>            writes prefix.pub and prefix.priv (bits 40 to 64)\n"
         "  enc <in> <out> <key.pub>          compress and encrypt a file\n"
         "  dec <in> <out> <key.priv>         decrypt and restore a file\n"
         "  attack [minbits maxbits]          factor toy keys with Pollard's rho\n"
         "  bench [out.csv]                   timing runs, written as CSV\n"
         "  serve [port]                      web interface on 127.0.0.1 (default 8080)");
}

/* ---------- commands ---------- */
static int cmd_find(const Genome *g) {
    MatchFn fn[3] = {kmp_search, bm_search, rk_search};
    const char *nm[3] = {"KMP", "Boyer-Moore", "Rabin-Karp"};
    size_t cnt[3];
    uint64_t sig[3];
    printf("genome %s, %llu bases\n", g->name, LLU(g->n));
    for (int i = 0; i < 3; i++) {
        TargetVec tv = {0};
        double t0 = now_ms();
        find_targets(g, fn[i], &tv);
        double ms = now_ms() - t0;
        cnt[i] = tv.n;
        sig[i] = targetvec_signature(&tv);
        printf("  %-12s %llu targets  signature %016llx  %.2f ms\n", nm[i], LLU(tv.n), LLU(sig[i]), ms);
        targetvec_free(&tv);
    }
    int same = cnt[0] == cnt[1] && cnt[1] == cnt[2] && sig[0] == sig[1] && sig[1] == sig[2];
    printf("all three agree: %s\n", same ? "yes" : "NO");
    return same ? 0 : 1;
}

static int cmd_offtarget(const Genome *g, const char *guide, int k) {
    uint64_t code;
    if (guide_encode(guide, &code) != 0) return fail("the guide must be exactly 20 letters from A, C, G, T.");
    OtIndex ix;
    double t0 = now_ms();
    ot_build(g, k, &ix);
    double build = now_ms() - t0;
    HitVec a = {0}, b = {0};
    t0 = now_ms();
    ot_query(&ix, code, &b);
    double tq = now_ms() - t0;
    t0 = now_ms();
    ot_brute(g, code, k, &a);
    double tb = now_ms() - t0;
    hitvec_sort(&a);
    hitvec_sort(&b);
    int same = a.n == b.n && hitvec_signature(&a) == hitvec_signature(&b);
    printf("%llu sites within %d mismatches (index %.3f ms, brute force %.2f ms, index build %.1f ms)\n",
           LLU(b.n), k, tq, tb, build);
    printf("index and brute force agree: %s\n", same ? "yes" : "NO");
    for (size_t i = 0; i < b.n && i < 40; i++) {
        char m[GUIDE_LEN + 1];
        for (int j = 0; j < GUIDE_LEN; j++) m[j] = (b.a[i].mask >> j) & 1 ? '^' : '.';
        m[GUIDE_LEN] = 0;
        printf("  %c %-10llu mm=%d  %s\n", b.a[i].strand, LLU(b.a[i].pos), b.a[i].mm, m);
    }
    hitvec_free(&a); hitvec_free(&b); ot_free(&ix);
    return same ? 0 : 1;
}

static void print_ranked(const char *tag, const Ranked *r) {
    printf("  %-6s %c %-9llu %s %s  score %6.2f  spec %6.2f  gc %2.0f%%  off-targets %d\n", tag, r->strand,
           LLU(r->pos), r->guide, r->pam, r->score, r->spec, r->gc * 100, r->nhits);
}

static int cmd_rank(const Genome *g, size_t from, size_t to, size_t K, int k) {
    OtIndex ix;
    ot_build(g, k, &ix);
    RankResult rr;
    if (rank_region(g, &ix, from, to, K, &rr) != 0) { ot_free(&ix); return fail("the region is too short or K is zero."); }
    printf("%llu candidates, %llu scored, %llu rejected (TTTT), median specificity %.2f\n",
           LLU(rr.ncand), LLU(rr.nscored), LLU(rr.nrej_tttt), rr.median_spec);
    for (size_t i = 0; i < rr.ntop; i++) { char t[32]; snprintf(t, sizeof t, "#%llu", LLU(i + 1)); print_ranked(t, &rr.top[i]); }
    for (size_t i = 0; i < rr.nworst; i++) print_ranked("worst", &rr.worst[i]);
    rank_free(&rr);
    ot_free(&ix);
    return 0;
}

static int cmd_cpg(const Genome *g) {
    IslandVec iv = {0};
    double t0 = now_ms();
    cpg_islands(g->seq, g->n, &iv);
    printf("%llu islands in %.2f ms\n", LLU(iv.n), now_ms() - t0);
    for (size_t i = 0; i < iv.n; i++)
        printf("  %llu..%llu  length %llu  GC %.0f%%  obs/exp CpG %.2f\n", LLU(iv.a[i].start), LLU(iv.a[i].end),
               LLU(iv.a[i].end - iv.a[i].start), iv.a[i].gc * 100, iv.a[i].oe);
    islandvec_free(&iv);
    return 0;
}

static char *make_seq(const char *arg, uint64_t seed) {
    static const char L[4] = {'A', 'C', 'G', 'T'};
    if (!strncmp(arg, "rand:", 5)) {
        unsigned long long n;
        if (parse_ull(arg + 5, &n) != 0 || n == 0 || n > 8000) return NULL;
        Rng r;
        rng_seed(&r, seed);
        char *s = xmalloc((size_t)n + 1);
        for (size_t i = 0; i < n; i++) s[i] = L[rng_below(&r, 4)];
        s[n] = 0;
        return s;
    }
    size_t n = strlen(arg);
    char *s = xmalloc(n + 1);
    for (size_t i = 0; i < n; i++) {
        if (base2(arg[i]) < 0) { free(s); return NULL; }
        s[i] = (char)(arg[i] & ~32);
    }
    s[n] = 0;
    return n ? s : (free(s), NULL);
}

static void print_alignment(const Alignment *al) {
    for (size_t i = 0; i < al->len; i += 60) {
        size_t m = al->len - i < 60 ? al->len - i : 60;
        printf("  %.*s\n  ", (int)m, al->a + i);
        for (size_t j = 0; j < m; j++) putchar(al->a[i + j] == al->b[i + j] ? '|' : ' ');
        printf("\n  %.*s\n\n", (int)m, al->b + i);
    }
}

static int cmd_align(const char *x, const char *y) {
    char *a = make_seq(x, 1), *b = make_seq(y, 2);
    if (!a || !b) { free(a); free(b); return fail("sequences must be letters A, C, G, T or rand:LENGTH (at most 8000)."); }
    size_t m = strlen(a), n = strlen(b);
    Alignment nw, hb;
    int rc_nw = nw_align(a, m, b, n, &nw);
    hirschberg_align(a, m, b, n, &hb);
    size_t swb;
    int sw = sw_score(a, m, b, n, &swb);
    if (rc_nw == 0) printf("Needleman-Wunsch score %d, table memory %llu bytes\n", nw.score, LLU(nw.bytes));
    else printf("Needleman-Wunsch skipped: sequences too long for the full table\n");
    printf("Hirschberg score %d, memory %llu bytes\n", hb.score, LLU(hb.bytes));
    printf("Smith-Waterman local score %d, memory %llu bytes\n", sw, LLU(swb));
    if (m + n <= 400) print_alignment(&hb);
    free(a); free(b);
    if (rc_nw == 0) alignment_free(&nw);
    alignment_free(&hb);
    return 0;
}

static int cmd_keygen(const char *bits, const char *prefix) {
    unsigned long long b;
    if (parse_ull(bits, &b) != 0 || b < 40 || b > 64) return fail("bits must be a number from 40 to 64.");
    Rng r;
    rng_seed(&r, os_random64());
    RsaKey k;
    rsa_keygen((int)b, &r, &k);
    char pub[512], priv[512], err[256];
    snprintf(pub, sizeof pub, "%s.pub", prefix);
    snprintf(priv, sizeof priv, "%s.priv", prefix);
    if (rsa_key_save(pub, &k, 0, err, sizeof err) || rsa_key_save(priv, &k, 1, err, sizeof err)) return fail(err);
    printf("wrote %s and %s (modulus %llu, e %llu)\n", pub, priv, LLU(k.n), LLU(k.e));
    return 0;
}

static int cmd_enc(const char *in, const char *out, const char *keyf) {
    char err[256];
    RsaKey k;
    if (rsa_key_load(keyf, &k, err, sizeof err) || k.e == 0) return fail(k.e == 0 && !*err ? "encryption needs the .pub key file." : err);
    Rng r;
    rng_seed(&r, os_random64());
    VaultStats st;
    if (vault_encrypt_file(in, out, &k, &r, &st, err, sizeof err)) return fail(err);
    printf("%llu bytes -> %llu compressed -> %llu encrypted\n", LLU(st.orig), LLU(st.comp), LLU(st.enc));
    return 0;
}

static int cmd_dec(const char *in, const char *out, const char *keyf) {
    char err[256];
    RsaKey k;
    if (rsa_key_load(keyf, &k, err, sizeof err)) return fail(err);
    if (k.d == 0) return fail("decryption needs the .priv key file.");
    if (vault_decrypt_file(in, out, &k, err, sizeof err)) return fail(err);
    puts("decrypted and checksum verified");
    return 0;
}

static void print_attack(const AttackRow *rows, int n) {
    puts("  prime bits  modulus bits  avg steps   theory 2^(b/2)  avg ms   recovered");
    for (int i = 0; i < n; i++)
        printf("  %-10d  %-12d  %-10.0f  %-14.0f  %-7.3f  %d of %d\n", rows[i].prime_bits, rows[i].prime_bits * 2,
               rows[i].avg_steps, rows[i].theory, rows[i].avg_ms, rows[i].ok, rows[i].trials);
}

static int cmd_attack(int minb, int maxb) {
    AttackRow rows[40];
    int n = attack_table(minb, maxb, 2, 20, 2024, rows, 40);
    print_attack(rows, n);
    if (n >= 2)
        printf("steps grew %.1fx from %d to %d bit primes; sqrt growth predicts %.1fx\n",
               rows[n - 1].avg_steps / rows[0].avg_steps, rows[0].prime_bits, rows[n - 1].prime_bits,
               (double)(1ULL << (rows[n - 1].prime_bits / 2)) / (double)(1ULL << (rows[0].prime_bits / 2)));
    return 0;
}

static int cmd_bench(const char *csv) {
    size_t sizes[4] = {100000, 300000, 1000000, 4600000};
    size_t lens[4] = {500, 1000, 2000, 4000};
    GenomeBenchRow gr[4];
    AlignBenchRow ar[4];
    if (bench_genomes(sizes, 4, gr) || bench_align(lens, 4, ar)) return fail("benchmark failed.");
    puts("  genome size   KMP ms    BM ms    RK ms   index build  index query  brute query  CpG ms");
    for (int i = 0; i < 4; i++)
        printf("  %-12llu  %-8.2f  %-7.2f  %-7.2f  %-11.1f  %-11.4f  %-11.2f  %.2f\n", LLU(gr[i].n), gr[i].kmp_ms,
               gr[i].bm_ms, gr[i].rk_ms, gr[i].build_ms, gr[i].query_ms, gr[i].brute_ms, gr[i].cpg_ms);
    puts("  align length  NW ms      Hirschberg ms   NW bytes      Hirschberg bytes");
    for (int i = 0; i < 4; i++)
        printf("  %-12llu  %-9.1f  %-14.1f  %-12llu  %llu\n", LLU(ar[i].len), ar[i].nw_ms, ar[i].hb_ms,
               LLU(ar[i].nw_bytes), LLU(ar[i].hb_bytes));
    if (csv) {
        FILE *f = fopen(csv, "wb");
        if (!f) return fail("cannot write the CSV file.");
        fputs("kind,size,kmp_ms,bm_ms,rk_ms,build_ms,query_ms,brute_ms,cpg_ms,sites\n", f);
        for (int i = 0; i < 4; i++)
            fprintf(f, "genome,%llu,%.3f,%.3f,%.3f,%.3f,%.5f,%.3f,%.3f,%llu\n", LLU(gr[i].n), gr[i].kmp_ms, gr[i].bm_ms,
                    gr[i].rk_ms, gr[i].build_ms, gr[i].query_ms, gr[i].brute_ms, gr[i].cpg_ms, LLU(gr[i].sites));
        fputs("kind,length,nw_ms,hirschberg_ms,nw_bytes,hirschberg_bytes\n", f);
        for (int i = 0; i < 4; i++)
            fprintf(f, "align,%llu,%.3f,%.3f,%llu,%llu\n", LLU(ar[i].len), ar[i].nw_ms, ar[i].hb_ms,
                    LLU(ar[i].nw_bytes), LLU(ar[i].hb_bytes));
        fclose(f);
        printf("wrote %s\n", csv);
    }
    return 0;
}

/* ---------- demo ---------- */
static int cmd_demo(size_t n) {
    Genome g;
    SynthInfo si;
    if (synth_genome(n, 42, &g, &si) != 0) return fail("demo genome needs at least 20000 bases.");
    printf("== Genome ==\n  %s, %llu bases. Planted guide %s at %llu with 12 decoys.\n", g.name, LLU(g.n), si.guide,
           LLU(si.guide_pos));

    puts("\n== 1. Target finder: every 20 letters + NGG on both strands ==");
    MatchFn fn[3] = {kmp_search, bm_search, rk_search};
    const char *nm[3] = {"KMP", "Boyer-Moore", "Rabin-Karp"};
    size_t cnt[3];
    uint64_t sig[3];
    for (int i = 0; i < 3; i++) {
        TargetVec tv = {0};
        double t0 = now_ms();
        find_targets(&g, fn[i], &tv);
        double ms = now_ms() - t0;
        cnt[i] = tv.n; sig[i] = targetvec_signature(&tv);
        printf("  %-12s %llu targets, %.2f ms\n", nm[i], LLU(tv.n), ms);
        targetvec_free(&tv);
    }
    check(cnt[0] == cnt[1] && cnt[1] == cnt[2] && sig[0] == sig[1] && sig[1] == sig[2], "three matchers: same count and signature");

    puts("\n== 2. Off-target scan: seed index against brute force (k = 3) ==");
    OtIndex ix;
    double t0 = now_ms();
    ot_build(&g, 3, &ix);
    double build = now_ms() - t0;
    uint64_t code;
    guide_encode(si.guide, &code);
    HitVec a = {0}, b = {0};
    t0 = now_ms(); ot_brute(&g, code, 3, &a); double tb = now_ms() - t0;
    t0 = now_ms(); ot_query(&ix, code, &b); double tq = now_ms() - t0;
    hitvec_sort(&a); hitvec_sort(&b);
    int by[4] = {0};
    for (size_t i = 0; i < b.n; i++) by[b.a[i].mm]++;
    printf("  planted guide: %llu sites (%d intended, %d with 1 mismatch, %d with 2, %d with 3)\n", LLU(b.n), by[0], by[1], by[2], by[3]);
    printf("  brute force %.3f ms, index query %.4f ms (build %.1f ms once)\n", tb, tq, build);
    check(a.n == b.n && hitvec_signature(&a) == hitvec_signature(&b), "index equals brute force");
    check(b.n == 13 && by[0] == 1, "1 intended site + 12 decoys");
    hitvec_free(&a); hitvec_free(&b);

    puts("\n== 3. Ranking: 30 kb around the planted guide ==");
    size_t from = si.guide_pos > 15000 ? si.guide_pos - 15000 : 0;
    RankResult rr;
    rank_region(&g, &ix, from, from + 30000, 5, &rr);
    printf("  %llu candidates, median specificity %.2f\n", LLU(rr.ncand), rr.median_spec);
    for (size_t i = 0; i < rr.ntop && i < 3; i++) print_ranked("best", &rr.top[i]);
    for (size_t i = 0; i < rr.nworst; i++) print_ranked("worst", &rr.worst[i]);
    int pw = 0;
    for (size_t i = 0; i < rr.nworst; i++) if (rr.worst[i].pos == si.guide_pos) pw = 1;
    check(pw, "planted guide is among the 3 worst");
    check(rr.ntop > 0 && rr.top[0].score == 100.0, "best guide has no near-matches (score 100)");
    rank_free(&rr);
    ot_free(&ix);

    puts("\n== 4. Epigenetics ==");
    IslandVec iv = {0};
    cpg_islands(g.seq, g.n, &iv);
    int hit = 0;
    for (int p = 0; p < 3; p++)
        for (size_t i = 0; i < iv.n; i++)
            if (iv.a[i].start < si.island_pos[p] + si.island_len[p] && iv.a[i].end > si.island_pos[p]) { hit++; break; }
    printf("  %llu islands found, %d of 3 planted islands overlapped\n", LLU(iv.n), hit);
    check(hit == 3, "all three planted CpG islands detected");
    islandvec_free(&iv);
    Genome h = {0};
    h.n = g.n; h.seq = xmalloc(g.n + 1); memcpy(h.seq, g.seq, g.n + 1);
    Rng rg; rng_seed(&rg, 3);
    genome_mutate(&h, 0.02, &rg);
    size_t L = 3000;
    size_t lcs = lcs_length(g.seq + 20000, L, h.seq + 20000, L);
    printf("  LCS of %llu letters vs the same region mutated at 2%%: %llu (%.1f%%)\n", LLU(L), LLU(lcs), 100.0 * (double)lcs / (double)L);
    check(lcs >= (size_t)(0.95 * (double)L) && lcs <= L, "LCS is close to the length for a 2% mutated copy");
    genome_free(&h);

    puts("\n== 5. Alignment (match +2, mismatch -1, gap -2) ==");
    size_t AL = 2000;
    char *x = xmalloc(AL + 1), *y = xmalloc(AL + 1);
    memcpy(x, g.seq + 50000, AL); x[AL] = 0;
    memcpy(y, x, AL + 1);
    Genome ty = {0}; ty.seq = y; ty.n = AL;
    rng_seed(&rg, 4);
    genome_mutate(&ty, 0.08, &rg);
    Alignment nw, hb;
    nw_align(x, AL, y, AL, &nw);
    hirschberg_align(x, AL, y, AL, &hb);
    size_t swb;
    int sw = sw_score(x, AL, y, AL, &swb);
    printf("  NW %d (%llu bytes), Hirschberg %d (%llu bytes), Smith-Waterman local %d\n", nw.score, LLU(nw.bytes), hb.score, LLU(hb.bytes), sw);
    check(nw.score == hb.score && alignment_score(&nw) == nw.score, "Hirschberg score equals full-table score");
    alignment_free(&nw); alignment_free(&hb);
    free(x); free(y);

    puts("\n== 6. Genome vault ==");
    Genome part = {0};
    part.n = g.n < 100000 ? g.n : 100000; part.seq = g.seq;
    snprintf(part.name, sizeof part.name, "demo_part");
    const char *f0 = "crispr_demo_plain.fa", *f1 = "crispr_demo_enc.bin", *f2 = "crispr_demo_back.fa";
    char err[256];
    Rng kr; rng_seed(&kr, 11);
    RsaKey k1, k2;
    rsa_keygen(64, &kr, &k1);
    rsa_keygen(64, &kr, &k2);
    genome_save(f0, &part, err, sizeof err);
    VaultStats st;
    int ok = vault_encrypt_file(f0, f1, &k1, &kr, &st, err, sizeof err) == 0;
    printf("  %llu bytes -> %llu compressed -> %llu encrypted\n", LLU(st.orig), LLU(st.comp), LLU(st.enc));
    ok = ok && vault_decrypt_file(f1, f2, &k1, err, sizeof err) == 0;
    check(ok && files_equal(f0, f2), "round trip: restored file is byte-identical");
    RsaKey bad = k2;
    check(vault_decrypt_file(f1, "crispr_demo_bad.fa", &bad, err, sizeof err) != 0, "different key pair rejected");
    bad = k1; bad.d = k1.d ^ 0x10;
    check(vault_decrypt_file(f1, "crispr_demo_bad.fa", &bad, err, sizeof err) != 0, "same modulus, wrong private exponent rejected");
    remove(f0); remove(f1); remove(f2); remove("crispr_demo_bad.fa");

    puts("\n== 7. Attack: Pollard's rho on toy keys ==");
    AttackRow rows[10];
    int nr = attack_table(12, 32, 4, 20, 2024, rows, 10);
    print_attack(rows, nr);
    int allok = 1;
    for (int i = 0; i < nr; i++) if (rows[i].ok != rows[i].trials) allok = 0;
    check(allok, "private key recovered for every 12 to 32 bit prime");

    genome_free(&g);
    printf("\n%d of %d checks passed.\n", checks - check_fail, checks);
    return check_fail ? 1 : 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { usage(); return 0; }
    const char *c = argv[1];
    Genome g = {0};
    unsigned long long v1, v2, v3;
    int rc;
    if (!strcmp(c, "demo")) {
        v1 = 300000;
        if (argc > 2 && parse_ull(argv[2], &v1) != 0) return fail("n must be a number.");
        return cmd_demo((size_t)v1);
    }
    if (!strcmp(c, "gen")) {
        if (argc < 4 || parse_ull(argv[2], &v1) != 0) return fail("usage: crispr gen <n> <out.fa> [seed]");
        v2 = 42;
        if (argc > 4 && parse_ull(argv[4], &v2) != 0) return fail("seed must be a number.");
        char err[256];
        if (synth_genome((size_t)v1, v2, &g, NULL) != 0) return fail("n must be at least 20000.");
        rc = genome_save(argv[3], &g, err, sizeof err);
        if (rc) fprintf(stderr, "error: %s\n", err); else printf("wrote %s (%llu bases)\n", argv[3], LLU(g.n));
        genome_free(&g);
        return rc != 0;
    }
    if (!strcmp(c, "find") || !strcmp(c, "cpg")) {
        if (argc < 3) return fail("a genome argument is needed.");
        if (load_genome_arg(argv[2], &g)) return 1;
        rc = !strcmp(c, "find") ? cmd_find(&g) : cmd_cpg(&g);
        genome_free(&g);
        return rc;
    }
    if (!strcmp(c, "offtarget")) {
        if (argc < 4) return fail("usage: crispr offtarget <genome> <guide> [k]");
        v1 = 3;
        if (argc > 4 && (parse_ull(argv[4], &v1) != 0 || v1 > 5)) return fail("k must be a number from 0 to 5.");
        if (load_genome_arg(argv[2], &g)) return 1;
        rc = cmd_offtarget(&g, argv[3], (int)v1);
        genome_free(&g);
        return rc;
    }
    if (!strcmp(c, "rank")) {
        if (argc < 5 || parse_ull(argv[3], &v1) || parse_ull(argv[4], &v2)) return fail("usage: crispr rank <genome> <from> <to> [K] [k]");
        v3 = 10;
        unsigned long long kk = 3;
        if (argc > 5 && (parse_ull(argv[5], &v3) || v3 == 0)) return fail("K must be a positive number.");
        if (argc > 6 && (parse_ull(argv[6], &kk) || kk > 5)) return fail("k must be a number from 0 to 5.");
        if (load_genome_arg(argv[2], &g)) return 1;
        rc = cmd_rank(&g, (size_t)v1, (size_t)v2, (size_t)v3, (int)kk);
        genome_free(&g);
        return rc;
    }
    if (!strcmp(c, "lcs")) {
        if (argc < 4) return fail("usage: crispr lcs <genomeA> <genomeB> [len]");
        Genome h = {0};
        v1 = 5000;
        if (argc > 4 && (parse_ull(argv[4], &v1) || v1 == 0)) return fail("len must be a positive number.");
        if (load_genome_arg(argv[2], &g)) return 1;
        if (load_genome_arg(argv[3], &h)) { genome_free(&g); return 1; }
        size_t la = g.n < v1 ? g.n : (size_t)v1, lb = h.n < v1 ? h.n : (size_t)v1;
        printf("LCS of the first %llu and %llu letters: %llu\n", LLU(la), LLU(lb), LLU(lcs_length(g.seq, la, h.seq, lb)));
        genome_free(&g); genome_free(&h);
        return 0;
    }
    if (!strcmp(c, "align")) { if (argc < 4) return fail("usage: crispr align <A> <B>"); return cmd_align(argv[2], argv[3]); }
    if (!strcmp(c, "keygen")) { if (argc < 4) return fail("usage: crispr keygen <bits> <prefix>"); return cmd_keygen(argv[2], argv[3]); }
    if (!strcmp(c, "enc")) { if (argc < 5) return fail("usage: crispr enc <in> <out> <key.pub>"); return cmd_enc(argv[2], argv[3], argv[4]); }
    if (!strcmp(c, "dec")) { if (argc < 5) return fail("usage: crispr dec <in> <out> <key.priv>"); return cmd_dec(argv[2], argv[3], argv[4]); }
    if (!strcmp(c, "attack")) {
        v1 = 12; v2 = 32;
        if (argc > 3 && (parse_ull(argv[2], &v1) || parse_ull(argv[3], &v2) || v1 < 8 || v2 > 32 || v1 > v2))
            return fail("bits must satisfy 8 <= minbits <= maxbits <= 32.");
        return cmd_attack((int)v1, (int)v2);
    }
    if (!strcmp(c, "bench")) return cmd_bench(argc > 2 ? argv[2] : NULL);
    if (!strcmp(c, "serve")) {
        v1 = 8080;
        if (argc > 2 && (parse_ull(argv[2], &v1) || v1 < 1 || v1 > 65535)) return fail("port must be 1 to 65535.");
        return serve_main((int)v1, "web");
    }
    fprintf(stderr, "error: unknown command '%s'. Run crispr with no arguments for the list.\n", c);
    return 1;
}
