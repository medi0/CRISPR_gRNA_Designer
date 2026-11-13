/* genome.c: FASTA parse/load/save, reverse complement, RNG, synthetic genome,
 * mutation, dynamic arrays, timer.
 *
 * Complexity: parse, save, mutate and the synthetic generator are O(n) in the
 * genome length. Dynamic array push is amortised O(1) (capacity doubles).
 */
#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif
#include "crispr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

/* ---------- helpers ---------- */
void *xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) { fputs("out of memory\n", stderr); exit(1); }
    return p;
}
void *xcalloc(size_t n, size_t sz) {
    void *p = calloc(n ? n : 1, sz ? sz : 1);
    if (!p) { fputs("out of memory\n", stderr); exit(1); }
    return p;
}
static void *grow(void *p, size_t *cap, size_t elem) {
    size_t nc = *cap ? *cap * 2 : 64;
    void *q = realloc(p, nc * elem);
    if (!q) { fputs("out of memory\n", stderr); exit(1); }
    *cap = nc;
    return q;
}

double now_ms(void) {
#ifdef _WIN32
    LARGE_INTEGER f, c;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&c);
    return 1000.0 * (double)c.QuadPart / (double)f.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
#endif
}

/* ---------- random numbers: splitmix64 ---------- */
#define GAMMA 0x9E3779B97F4A7C15ULL
static uint64_t fmix(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
uint64_t mix64(uint64_t x) { return fmix(x + GAMMA); }
void rng_seed(Rng *r, uint64_t seed) { r->state = seed; }
uint64_t rng_next(Rng *r) { r->state += GAMMA; return fmix(r->state); }
uint32_t rng_below(Rng *r, uint32_t n) {
    return (uint32_t)(((rng_next(r) >> 32) * (uint64_t)n) >> 32);
}
double rng_unit(Rng *r) { return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0); }

/* ---------- dynamic arrays ---------- */
void posvec_push(PosVec *v, size_t x) {
    if (v->n == v->cap) v->a = grow(v->a, &v->cap, sizeof *v->a);
    v->a[v->n++] = x;
}
void posvec_free(PosVec *v) { free(v->a); v->a = NULL; v->n = v->cap = 0; }

void targetvec_push(TargetVec *v, Target t) {
    if (v->n == v->cap) v->a = grow(v->a, &v->cap, sizeof *v->a);
    v->a[v->n++] = t;
}
void targetvec_free(TargetVec *v) { free(v->a); v->a = NULL; v->n = v->cap = 0; }
uint64_t targetvec_signature(const TargetVec *v) {
    uint64_t s = 0;
    for (size_t i = 0; i < v->n; i++)
        s += mix64((uint64_t)v->a[i].pos * 2 + (v->a[i].strand == '-'));
    return s;
}

void hitvec_push(HitVec *v, Hit h) {
    if (v->n == v->cap) v->a = grow(v->a, &v->cap, sizeof *v->a);
    v->a[v->n++] = h;
}
void hitvec_free(HitVec *v) { free(v->a); v->a = NULL; v->n = v->cap = 0; }
static int hit_cmp(const void *x, const void *y) {
    const Hit *a = x, *b = y;
    if (a->pos != b->pos) return a->pos < b->pos ? -1 : 1;
    return (int)a->strand - (int)b->strand;
}
void hitvec_sort(HitVec *v) { if (v->n > 1) qsort(v->a, v->n, sizeof *v->a, hit_cmp); }
uint64_t hitvec_signature(const HitVec *v) {
    uint64_t s = 0;
    for (size_t i = 0; i < v->n; i++) {
        const Hit *h = &v->a[i];
        s += mix64((uint64_t)h->pos * 2 + (h->strand == '-'))
           ^ mix64(((uint64_t)h->mask << 4) | (uint64_t)h->mm);
    }
    return s;
}

/* ---------- sequence basics ---------- */
char complement(char c) {
    switch (c) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default: return 'N';
    }
}
void revcomp(const char *src, size_t n, char *dst) {
    for (size_t i = 0; i < n; i++) dst[i] = complement(src[n - 1 - i]);
}

/* ---------- FASTA ---------- */
static void set_err(char *err, size_t errsz, const char *msg) {
    if (err && errsz) snprintf(err, errsz, "%s", msg);
}

void genome_free(Genome *g) { free(g->seq); g->seq = NULL; g->n = 0; g->name[0] = 0; }

/* Reads the FIRST record only. Upper-cases, keeps only A C G T, which also
 * drops carriage returns, newlines, N and any other IUPAC letter. */
int genome_parse(const char *buf, size_t len, Genome *g, char *err, size_t errsz) {
    memset(g, 0, sizeof *g);
    size_t i = 0;
    while (i < len && (buf[i] == ' ' || buf[i] == '\t' || buf[i] == '\r' || buf[i] == '\n')) i++;
    if (i >= len || buf[i] != '>') {
        set_err(err, errsz, "This does not look like a FASTA file: the first line must start with '>'.");
        return -1;
    }
    i++;
    size_t k = 0;
    while (i < len && buf[i] != '\n' && buf[i] != '\r') {
        if (k < sizeof g->name - 1) g->name[k++] = buf[i];
        i++;
    }
    g->name[k] = 0;
    char *seq = xmalloc(len + 1);
    size_t n = 0;
    for (; i < len && buf[i] != '>'; i++) {
        char c = buf[i];
        if (c >= 'a' && c <= 'z') c = (char)(c - 32);
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') seq[n++] = c;
    }
    if (n == 0) {
        free(seq);
        set_err(err, errsz, "The FASTA record has no A, C, G or T letters.");
        return -1;
    }
    seq[n] = 0;
    g->seq = seq;
    g->n = n;
    return 0;
}

int genome_load(const char *path, Genome *g, char *err, size_t errsz) {
    FILE *f = fopen(path, "rb");
    if (!f) {
        if (err && errsz) snprintf(err, errsz, "Cannot open file '%s'.", path);
        return -1;
    }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz <= 0) { fclose(f); set_err(err, errsz, "The file is empty."); return -1; }
    char *buf = xmalloc((size_t)sz);
    size_t got = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    int rc = genome_parse(buf, got, g, err, errsz);
    free(buf);
    return rc;
}

int genome_save(const char *path, const Genome *g, char *err, size_t errsz) {
    FILE *f = fopen(path, "wb");
    if (!f) {
        if (err && errsz) snprintf(err, errsz, "Cannot write file '%s'.", path);
        return -1;
    }
    fprintf(f, ">%s\n", g->name[0] ? g->name : "genome");
    for (size_t i = 0; i < g->n; i += 70) {
        size_t m = g->n - i < 70 ? g->n - i : 70;
        fwrite(g->seq + i, 1, m, f);
        fputc('\n', f);
    }
    fclose(f);
    return 0;
}

/* ---------- synthetic genome ---------- */
static char draw(Rng *r, const double cum[3]) {
    double u = rng_unit(r);
    if (u < cum[0]) return 'A';
    if (u < cum[1]) return 'C';
    if (u < cum[2]) return 'G';
    return 'T';
}

static int overlaps(const size_t *rs, const size_t *re, int cnt, size_t a, size_t b) {
    for (int i = 0; i < cnt; i++)
        if (a < re[i] && b > rs[i]) return 1;
    return 0;
}

static void plant_site(char *seq, size_t p, char strand, const char *proto, char pam_n) {
    char site[SITE_LEN];
    memcpy(site, proto, GUIDE_LEN);
    site[20] = pam_n; site[21] = 'G'; site[22] = 'G';
    if (strand == '+') memcpy(seq + p, site, SITE_LEN);
    else revcomp(site, SITE_LEN, seq + p);
}

/* Background: AT-rich (about 36% GC) with CpG depleted (85% of CG pairs are
 * redrawn), like real genomes. Three CpG islands of 1200 letters. One planted
 * guide plus 12 decoys that differ from it by 1 to 3 letters and have a valid
 * NGG PAM, alternating strands. */
int synth_genome(size_t n, uint64_t seed, Genome *g, SynthInfo *info) {
    static const double BG[3] = {0.32, 0.50, 0.68};
    static const double ISL[3] = {0.18, 0.50, 0.82};
    static const char GUIDE[] = "GAGTCCGAGCAGAAGAAGAA";
    static const char LET[4] = {'A', 'C', 'G', 'T'};
    const size_t ISLAND_LEN = 1200;
    SynthInfo tmp;
    SynthInfo *si = info ? info : &tmp;
    memset(g, 0, sizeof *g);
    memset(si, 0, sizeof *si);
    if (n < 20000 || n > (size_t)0xF0000000u) return -1;

    Rng r;
    rng_seed(&r, seed);
    char *seq = xmalloc(n + 1);
    char prev = 'A';
    for (size_t i = 0; i < n; i++) {
        char b;
        do { b = draw(&r, BG); } while (prev == 'C' && b == 'G' && rng_unit(&r) < 0.85);
        seq[i] = b;
        prev = b;
    }
    seq[n] = 0;

    size_t rs[16], re[16];
    int nres = 0;
    size_t ipos[3] = {n / 5, n / 2, n / 5 * 4};
    for (int i = 0; i < 3; i++) {
        for (size_t j = 0; j < ISLAND_LEN; j++) seq[ipos[i] + j] = draw(&r, ISL);
        si->island_pos[i] = ipos[i];
        si->island_len[i] = ISLAND_LEN;
        rs[nres] = ipos[i]; re[nres] = ipos[i] + ISLAND_LEN; nres++;
    }

    size_t p0 = n / 3;
    memcpy(si->guide, GUIDE, GUIDE_LEN + 1);
    si->guide_pos = p0;
    plant_site(seq, p0, '+', GUIDE, 'C');
    rs[nres] = p0; re[nres] = p0 + SITE_LEN; nres++;

    for (int i = 0; i < 12; i++) {
        size_t p;
        do {
            p = 100 + rng_below(&r, (uint32_t)(n - 200 - SITE_LEN));
        } while (overlaps(rs, re, nres, p - 50, p + SITE_LEN + 50));
        char proto[GUIDE_LEN + 1];
        memcpy(proto, GUIDE, GUIDE_LEN + 1);
        int m = 1 + (i % 3);
        int used[GUIDE_LEN] = {0};
        for (int c = 0; c < m; c++) {
            int pos;
            do { pos = (int)rng_below(&r, GUIDE_LEN); } while (used[pos]);
            used[pos] = 1;
            int orig = base2(proto[pos]);
            proto[pos] = LET[(orig + 1 + (int)rng_below(&r, 3)) % 4];
        }
        char strand = (i % 2 == 0) ? '+' : '-';
        plant_site(seq, p, strand, proto, LET[rng_below(&r, 4)]);
        rs[nres] = p; re[nres] = p + SITE_LEN; nres++;
        si->decoy_pos[i] = p;
        si->decoy_strand[i] = strand;
        si->decoy_mm[i] = m;
    }

    g->seq = seq;
    g->n = n;
    snprintf(g->name, sizeof g->name, "synthetic_%zu_seed%llu", n, (unsigned long long)seed);
    return 0;
}

void genome_mutate(Genome *g, double rate, Rng *r) {
    static const char LET[4] = {'A', 'C', 'G', 'T'};
    for (size_t i = 0; i < g->n; i++)
        if (rng_unit(r) < rate)
            g->seq[i] = LET[(base2(g->seq[i]) + 1 + (int)rng_below(r, 3)) % 4];
}
