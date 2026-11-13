/* crispr.h: shared types and prototypes for the CRISPR Guide RNA Designer.
 *
 * Conventions used everywhere:
 *  - Coordinates are 0-based. A "site" is a 23-letter window on the forward
 *    strand starting at `pos`: 20-letter protospacer + 3-letter PAM.
 *  - strand '+': protospacer = seq[pos .. pos+19], PAM = seq[pos+20 .. pos+22],
 *    and the PAM must look like NGG, so seq[pos+21] = seq[pos+22] = 'G'.
 *  - strand '-': the window reads CCN on the forward strand, so
 *    seq[pos] = seq[pos+1] = 'C'. The guide is the reverse complement of
 *    seq[pos+3 .. pos+22].
 *  - Guide letter index 0 is PAM-distal, index 19 touches the PAM.
 *  - Modules return data. Only main.c prints and only server.c builds JSON.
 */
#ifndef CRISPR_H
#define CRISPR_H

#include <stddef.h>
#include <stdint.h>

#define GUIDE_LEN 20
#define SITE_LEN 23

/* ---------- small helpers ---------- */
void *xmalloc(size_t n);
void *xcalloc(size_t n, size_t sz);
double now_ms(void);
uint64_t mix64(uint64_t x);

static inline int base2(char c) { /* A0 C1 G2 T3 (complement = 3 - code) */
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
    }
}

/* ---------- random numbers (splitmix64, deterministic) ---------- */
typedef struct { uint64_t state; } Rng;
void rng_seed(Rng *r, uint64_t seed);
uint64_t rng_next(Rng *r);
uint32_t rng_below(Rng *r, uint32_t n);
double rng_unit(Rng *r);

/* ---------- dynamic arrays ---------- */
typedef struct { size_t *a; size_t n, cap; } PosVec;
void posvec_push(PosVec *v, size_t x);
void posvec_free(PosVec *v);

typedef struct { size_t pos; char strand; } Target;
typedef struct { Target *a; size_t n, cap; } TargetVec;
void targetvec_push(TargetVec *v, Target t);
void targetvec_free(TargetVec *v);
uint64_t targetvec_signature(const TargetVec *v); /* order independent */

typedef struct { size_t pos; char strand; int mm; uint32_t mask; } Hit;
typedef struct { Hit *a; size_t n, cap; } HitVec;
void hitvec_push(HitVec *v, Hit h);
void hitvec_free(HitVec *v);
void hitvec_sort(HitVec *v);                  /* by (pos, strand) */
uint64_t hitvec_signature(const HitVec *v);   /* order independent */

/* ---------- genome.c ---------- */
typedef struct { char *seq; size_t n; char name[128]; } Genome;

typedef struct {
    char guide[GUIDE_LEN + 1];   /* planted guide letters */
    size_t guide_pos;            /* its site (strand '+') */
    size_t decoy_pos[12];
    char decoy_strand[12];
    int decoy_mm[12];            /* 1 to 3 mismatches each */
    size_t island_pos[3];
    size_t island_len[3];
} SynthInfo;

char complement(char c);
void revcomp(const char *src, size_t n, char *dst); /* dst must not overlap src */
void genome_free(Genome *g);
int genome_parse(const char *buf, size_t len, Genome *g, char *err, size_t errsz);
int genome_load(const char *path, Genome *g, char *err, size_t errsz);
int genome_save(const char *path, const Genome *g, char *err, size_t errsz);
int synth_genome(size_t n, uint64_t seed, Genome *g, SynthInfo *info);
void genome_mutate(Genome *g, double rate, Rng *r);

/* ---------- strmatch.c ---------- */
/* Each matcher counts every (overlapping) occurrence; if out != NULL the
 * start positions are appended to it. */
typedef size_t (*MatchFn)(const char *t, size_t n, const char *p, size_t m, PosVec *out);
size_t kmp_search(const char *t, size_t n, const char *p, size_t m, PosVec *out);
size_t bm_search(const char *t, size_t n, const char *p, size_t m, PosVec *out);
size_t rk_search(const char *t, size_t n, const char *p, size_t m, PosVec *out);
void find_targets(const Genome *g, MatchFn fn, TargetVec *out);

/* ---------- offtarget.c ---------- */
typedef struct { size_t pos; char strand; uint64_t code; } Site;
typedef struct {
    int k, nseg;
    int seg_start[6], seg_len[6];
    size_t nsites, nbuckets;
    Site *sites;
    uint32_t *off[6];   /* CSR offsets, nbuckets + 1 entries */
    uint32_t *ent[6];   /* site indices, nsites entries */
} OtIndex;

int guide_encode(const char *s, uint64_t *code);   /* 0 ok, -1 bad guide */
void guide_decode(uint64_t code, char out[GUIDE_LEN + 1]);
uint64_t site_code(const char *seq, size_t p, char strand);
int ot_build(const Genome *g, int k, OtIndex *ix);
void ot_free(OtIndex *ix);
void ot_query(const OtIndex *ix, uint64_t guide, HitVec *out);
void ot_brute(const Genome *g, uint64_t guide, int k, HitVec *out);

/* ---------- scoring.c ---------- */
typedef struct {
    size_t pos; char strand;
    char guide[GUIDE_LEN + 1];
    char pam[4];
    double gc;      /* fraction 0..1 */
    double risk;    /* sum of per-hit risks, on-target excluded */
    double spec;    /* 100 / (1 + risk) */
    double score;   /* spec minus 20 if GC outside 40-70% */
    int nhits;      /* off-target sites within k mismatches */
} Ranked;

typedef struct {
    Ranked *top; size_t ntop;
    Ranked worst[3]; size_t nworst;
    double median_spec;
    size_t ncand, nscored, nrej_tttt;
} RankResult;

double hit_risk(const Hit *h);
int rank_region(const Genome *g, const OtIndex *ix, size_t from, size_t to,
                size_t K, RankResult *res);
void rank_free(RankResult *res);

/* ---------- epigen.c ---------- */
typedef struct { size_t start, end; double gc, oe; } CpgIsland; /* [start, end) */
typedef struct { CpgIsland *a; size_t n, cap; } IslandVec;
void islandvec_free(IslandVec *v);
void cpg_islands(const char *seq, size_t n, IslandVec *out);
void gc_track(const char *seq, size_t n, size_t bins, double *out); /* GC fraction per bin */
size_t lcs_length(const char *a, size_t la, const char *b, size_t lb);

/* ---------- align.c ---------- */
#define AL_MATCH 2
#define AL_MISMATCH (-1)
#define AL_GAP (-2)
typedef struct {
    int score;
    char *a, *b;     /* aligned strings, '-' marks a gap, NUL terminated */
    size_t len;
    size_t bytes;    /* working memory the algorithm used */
} Alignment;
int nw_align(const char *a, size_t m, const char *b, size_t n, Alignment *out);        /* -1 if m+n > 16000 */
int hirschberg_align(const char *a, size_t m, const char *b, size_t n, Alignment *out);
int sw_score(const char *a, size_t m, const char *b, size_t n, size_t *bytes);
int alignment_score(const Alignment *al);
void alignment_free(Alignment *al);
size_t peak_rss_bytes(void);

/* ---------- vault.c ---------- */
typedef struct { uint64_t n, e, d, p, q; int bits; } RsaKey;
typedef struct { size_t orig, comp, enc; } VaultStats;
typedef struct { int prime_bits, trials, ok; double avg_steps, theory, avg_ms; } AttackRow;
uint64_t os_random64(void);
uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m);
uint64_t powmod(uint64_t b, uint64_t e, uint64_t m);
int is_prime64(uint64_t n);
uint64_t fnv1a(const uint8_t *p, size_t n);
int rsa_keygen(int modbits, Rng *r, RsaKey *k);   /* modbits 16..64, 0 on success */
int rsa_key_save(const char *path, const RsaKey *k, int priv, char *err, size_t errsz);
int rsa_key_load(const char *path, RsaKey *k, char *err, size_t errsz);
int huff_compress(const uint8_t *in, size_t n, uint8_t **out, size_t *outn);
int huff_decompress(const uint8_t *in, size_t n, uint8_t **out, size_t *outn);
int vault_encrypt(const uint8_t *in, size_t n, const RsaKey *k, Rng *r,
                  uint8_t **out, size_t *outn, VaultStats *st);
int vault_decrypt(const uint8_t *in, size_t n, const RsaKey *k,
                  uint8_t **out, size_t *outn, char *err, size_t errsz);
int vault_encrypt_file(const char *in, const char *out, const RsaKey *k, Rng *r,
                       VaultStats *st, char *err, size_t errsz);
int vault_decrypt_file(const char *in, const char *out, const RsaKey *k, char *err, size_t errsz);
int files_equal(const char *p1, const char *p2);   /* 1 equal, 0 different or unreadable */
uint64_t pollard_rho(uint64_t n, uint64_t *steps);
int attack_table(int minb, int maxb, int step, int trials, uint64_t seed, AttackRow *rows, int cap);

/* ---------- bench.c ---------- */
typedef struct {
    size_t n, sites;
    double kmp_ms, bm_ms, rk_ms, build_ms, query_ms, brute_ms, cpg_ms;
} GenomeBenchRow;
typedef struct { size_t len; double nw_ms, hb_ms; size_t nw_bytes, hb_bytes; } AlignBenchRow;
int bench_genomes(const size_t *sizes, int ns, GenomeBenchRow *rows);
int bench_align(const size_t *lens, int nl, AlignBenchRow *rows);

/* ---------- server.c ---------- */
int serve_main(int port, const char *webdir);

#endif
