/* server.c: single-threaded blocking HTTP server on 127.0.0.1 and the JSON API.
 *
 * One genome and one cached off-target index live in memory. The index is
 * rebuilt only when k changes or a new genome is loaded. All JSON is built
 * here; the algorithms live in the other modules and never print.
 * Positions are 0-based and mark the start of the 23-letter site.
 */
#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif
#include "crispr.h"
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <winsock2.h>
typedef SOCKET sock_t;
#define CLOSESOCK closesocket
#define BADSOCK INVALID_SOCKET
#else
#include <arpa/inet.h>
#include <netinet/in.h>
#include <signal.h>
#include <sys/socket.h>
#include <unistd.h>
typedef int sock_t;
#define CLOSESOCK close
#define BADSOCK (-1)
#endif

#define LLU(x) ((unsigned long long)(x))
#define MAX_BODY (64u * 1024u * 1024u)

/* ---------- state ---------- */
static Genome G;
static int loaded = 0;
static SynthInfo SI;
static int has_si = 0;
static OtIndex IX;
static int ix_k = -1;

static void set_genome(Genome *g, const SynthInfo *si) {
    if (ix_k >= 0) { ot_free(&IX); ix_k = -1; }
    if (loaded) genome_free(&G);
    G = *g;
    loaded = 1;
    has_si = si != NULL;
    if (si) SI = *si;
}

static const OtIndex *get_index(int k, double *build_ms) {
    *build_ms = 0;
    if (ix_k != k) {
        if (ix_k >= 0) ot_free(&IX);
        double t0 = now_ms();
        ot_build(&G, k, &IX);
        *build_ms = now_ms() - t0;
        ix_k = k;
    }
    return &IX;
}

/* ---------- growable text buffer ---------- */
typedef struct { char *s; size_t n, cap; } Buf;

static void bputs(Buf *b, const char *s) {
    size_t l = strlen(s);
    if (b->n + l + 1 > b->cap) {
        b->cap = (b->n + l + 1) * 2;
        b->s = realloc(b->s, b->cap);
        if (!b->s) exit(1);
    }
    memcpy(b->s + b->n, s, l + 1);
    b->n += l;
}
/* GCC's -Wformat-truncation misfires on a wrapper whose format is a parameter
 * (seen with -O1 and sanitizers). Every call site passes a literal format. */
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-truncation"
#endif
static void bprintf(Buf *b, const char *fmt, ...) {
    char tmp[1024];
    va_list ap;
    va_start(ap, fmt);
    int l = vsnprintf(tmp, sizeof tmp, fmt, ap);
    va_end(ap);
    if (l < (int)sizeof tmp) { bputs(b, tmp); return; }
    char *big = xmalloc((size_t)l + 1);
    va_start(ap, fmt);
    vsnprintf(big, (size_t)l + 1, fmt, ap);
    va_end(ap);
    bputs(b, big);
    free(big);
}
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
static void bstr(Buf *b, const char *s) { /* JSON string with escapes */
    bputs(b, "\"");
    for (; *s; s++) {
        unsigned char c = (unsigned char)*s;
        if (c == '"' || c == '\\') { char t[3] = {'\\', (char)c, 0}; bputs(b, t); }
        else if (c < 32) { char t[8]; snprintf(t, sizeof t, "\\u%04x", c); bputs(b, t); }
        else { char t[2] = {(char)c, 0}; bputs(b, t); }
    }
    bputs(b, "\"");
}
static double fin(double x) { return isfinite(x) ? x : 0.0; }
static double ratio(double a, double b) { return b > 1e-9 ? a / b : 0.0; }

static int jerr(Buf *o, int *status, const char *msg) {
    o->n = 0;
    if (o->s) o->s[0] = 0;
    bputs(o, "{\"error\":");
    bstr(o, msg);
    bputs(o, "}");
    *status = 400;
    return 0;
}

/* ---------- query parameters ---------- */
static int qstr(const char *q, const char *key, char *out, size_t sz) {
    size_t kl = strlen(key);
    while (q && *q) {
        if (!strncmp(q, key, kl) && q[kl] == '=') {
            const char *v = q + kl + 1;
            size_t i = 0;
            while (*v && *v != '&') {
                char c = *v;
                if (c == '+') c = ' ';
                else if (c == '%' && v[1] && v[2]) {
                    char h[3] = {v[1], v[2], 0};
                    c = (char)strtol(h, NULL, 16);
                    v += 2;
                }
                if (i + 1 < sz) out[i++] = c;
                v++;
            }
            out[i] = 0;
            return 1;
        }
        q = strchr(q, '&');
        if (q) q++;
    }
    return 0;
}

/* 0 ok (value or default), -1 invalid */
static int qint(const char *q, const char *key, long long def, long long lo, long long hi, long long *v) {
    char t[64], *end;
    if (!qstr(q, key, t, sizeof t)) { *v = def; return 0; }
    if (!*t) return -1;
    *v = strtoll(t, &end, 10);
    return (*end || *v < lo || *v > hi) ? -1 : 0;
}

static int bad_param(Buf *o, int *st, const char *name, long long lo, long long hi) {
    char m[160];
    snprintf(m, sizeof m, "Parameter '%s' must be a whole number from %lld to %lld.", name, lo, hi);
    return jerr(o, st, m);
}

static int need_genome(Buf *o, int *st) {
    if (loaded) return 0;
    jerr(o, st, "No genome is loaded. Upload a FASTA file or pick a synthetic genome on the Genome page.");
    return 1;
}

/* ---------- endpoints ---------- */
static void status_json(Buf *o) {
    bputs(o, "{\"loaded\":");
    if (!loaded) { bputs(o, "false,\"index_k\":-1}"); return; }
    size_t gc = 0;
    for (size_t i = 0; i < G.n; i++) gc += (G.seq[i] == 'G' || G.seq[i] == 'C');
    bprintf(o, "true,\"name\":");
    bstr(o, G.name);
    bprintf(o, ",\"n\":%llu,\"gc\":%.4f,\"index_k\":%d,\"planted\":", LLU(G.n), (double)gc / (double)G.n, ix_k);
    if (has_si) {
        bprintf(o, "{\"guide\":\"%s\",\"pos\":%llu,\"decoys\":12,\"islands\":[", SI.guide, LLU(SI.guide_pos));
        for (int i = 0; i < 3; i++) bprintf(o, "%s{\"start\":%llu,\"len\":%llu}", i ? "," : "", LLU(SI.island_pos[i]), LLU(SI.island_len[i]));
        bputs(o, "]}");
    } else bputs(o, "null");
    bputs(o, "}");
}

static int api_synth(const char *q, Buf *o, int *st) {
    long long n;
    if (qint(q, "n", 100000, 20000, 10000000, &n)) return bad_param(o, st, "n", 20000, 10000000);
    Genome g;
    SynthInfo si;
    if (synth_genome((size_t)n, 42, &g, &si)) return jerr(o, st, "Could not build the synthetic genome.");
    set_genome(&g, &si);
    status_json(o);
    return 0;
}

static int api_upload(const char *q, const char *body, size_t blen, Buf *o, int *st) {
    char err[256], name[128];
    if (blen == 0) return jerr(o, st, "No file was received.");
    Genome g;
    if (genome_parse(body, blen, &g, err, sizeof err)) return jerr(o, st, err);
    if (qstr(q, "name", name, sizeof name) && *name) snprintf(g.name, sizeof g.name, "%s", name);
    set_genome(&g, NULL);
    status_json(o);
    return 0;
}

static int api_find(Buf *o, int *st) {
    if (need_genome(o, st)) return 0;
    MatchFn fn[3] = {kmp_search, bm_search, rk_search};
    const char *nm[3] = {"KMP", "Boyer-Moore", "Rabin-Karp"};
    size_t cnt[3];
    uint64_t sig[3];
    double ms[3];
    for (int i = 0; i < 3; i++) {
        TargetVec tv = {0};
        double t0 = now_ms();
        find_targets(&G, fn[i], &tv);
        ms[i] = now_ms() - t0;
        cnt[i] = tv.n;
        sig[i] = targetvec_signature(&tv);
        targetvec_free(&tv);
    }
    int agree = cnt[0] == cnt[1] && cnt[1] == cnt[2] && sig[0] == sig[1] && sig[1] == sig[2];
    bprintf(o, "{\"genome_n\":%llu,\"agree\":%s,\"matchers\":[", LLU(G.n), agree ? "true" : "false");
    for (int i = 0; i < 3; i++)
        bprintf(o, "%s{\"name\":\"%s\",\"targets\":%llu,\"signature\":\"%016llx\",\"ms\":%.3f}", i ? "," : "", nm[i], LLU(cnt[i]), LLU(sig[i]), fin(ms[i]));
    bputs(o, "],\"pattern16\":");
    if (G.n >= 1000) {
        size_t at = has_si ? SI.guide_pos : G.n / 3;
        char pat[17];
        memcpy(pat, G.seq + at, 16);
        pat[16] = 0;
        size_t pc[3];
        double pm[3];
        for (int i = 0; i < 3; i++) {
            double t0 = now_ms();
            for (int r = 0; r < 3; r++) pc[i] = fn[i](G.seq, G.n, pat, 16, NULL);
            pm[i] = (now_ms() - t0) / 3;
        }
        bprintf(o, "{\"pattern\":\"%s\",\"agree\":%s,\"results\":[", pat, (pc[0] == pc[1] && pc[1] == pc[2]) ? "true" : "false");
        for (int i = 0; i < 3; i++)
            bprintf(o, "%s{\"name\":\"%s\",\"matches\":%llu,\"ms\":%.3f}", i ? "," : "", nm[i], LLU(pc[i]), fin(pm[i]));
        bputs(o, "]}");
    } else bputs(o, "null");
    bputs(o, "}");
    return 0;
}

static void ranked_json(Buf *o, const Ranked *r) {
    bprintf(o, "{\"pos\":%llu,\"strand\":\"%c\",\"guide\":\"%s\",\"pam\":\"%s\",\"gc\":%.3f,\"spec\":%.2f,\"score\":%.2f,\"hits\":%d}",
            LLU(r->pos), r->strand, r->guide, r->pam, r->gc, fin(r->spec), fin(r->score), r->nhits);
}

static int api_rank(const char *q, Buf *o, int *st) {
    if (need_genome(o, st)) return 0;
    long long from, to, K, k;
    long long def_from = (has_si && SI.guide_pos > 15000) ? (long long)SI.guide_pos - 15000 : 0;
    if (qint(q, "from", def_from, 0, (long long)G.n, &from)) return bad_param(o, st, "from", 0, (long long)G.n);
    if (qint(q, "to", from + 30000 < (long long)G.n ? from + 30000 : (long long)G.n, 0, (long long)G.n, &to))
        return bad_param(o, st, "to", 0, (long long)G.n);
    if (qint(q, "K", 10, 1, 100, &K)) return bad_param(o, st, "K", 1, 100);
    if (qint(q, "k", 3, 0, 5, &k)) return bad_param(o, st, "k", 0, 5);
    if (to - from < SITE_LEN) return jerr(o, st, "The region must be at least 23 bases long: 'to' has to be larger than 'from'.");
    if (to - from > 30000) return jerr(o, st, "The region is limited to 30,000 bases. Pick a shorter stretch.");
    double bms;
    const OtIndex *ix = get_index((int)k, &bms);
    RankResult rr;
    double t0 = now_ms();
    if (rank_region(&G, ix, (size_t)from, (size_t)to, (size_t)K, &rr)) return jerr(o, st, "Could not rank this region.");
    double ms = now_ms() - t0;
    bprintf(o, "{\"from\":%lld,\"to\":%lld,\"K\":%lld,\"k\":%lld,\"candidates\":%llu,\"scored\":%llu,\"rejected_tttt\":%llu,"
               "\"median_spec\":%.2f,\"ms\":%.2f,\"index_ms\":%.2f,\"planted_pos\":",
            from, to, K, k, LLU(rr.ncand), LLU(rr.nscored), LLU(rr.nrej_tttt), fin(rr.median_spec), fin(ms), fin(bms));
    if (has_si) bprintf(o, "%llu", LLU(SI.guide_pos)); else bputs(o, "null");
    bputs(o, ",\"top\":[");
    for (size_t i = 0; i < rr.ntop; i++) { if (i) bputs(o, ","); ranked_json(o, &rr.top[i]); }
    bputs(o, "],\"worst\":[");
    for (size_t i = 0; i < rr.nworst; i++) { if (i) bputs(o, ","); ranked_json(o, &rr.worst[i]); }
    bputs(o, "]}");
    rank_free(&rr);
    return 0;
}

static int hit_order(const void *x, const void *y) {
    const Hit *a = x, *b = y;
    if (a->mm != b->mm) return a->mm - b->mm;
    if (a->pos != b->pos) return a->pos < b->pos ? -1 : 1;
    return (int)a->strand - (int)b->strand;
}

static int api_hits(const char *q, Buf *o, int *st) {
    if (need_genome(o, st)) return 0;
    char guide[64];
    uint64_t code;
    long long k;
    if (!qstr(q, "guide", guide, sizeof guide) || guide_encode(guide, &code))
        return jerr(o, st, "The guide must be exactly 20 letters from A, C, G, T.");
    if (qint(q, "k", 3, 0, 5, &k)) return bad_param(o, st, "k", 0, 5);
    double bms;
    const OtIndex *ix = get_index((int)k, &bms);
    HitVec hv = {0};
    double t0 = now_ms();
    ot_query(ix, code, &hv);
    double ms = now_ms() - t0;
    qsort(hv.a, hv.n, sizeof *hv.a, hit_order);
    char up[GUIDE_LEN + 1];
    guide_decode(code, up);
    bprintf(o, "{\"guide\":\"%s\",\"k\":%lld,\"count\":%llu,\"ms\":%.4f,\"index_ms\":%.2f,\"truncated\":%s,\"hits\":[", up, k,
            LLU(hv.n), fin(ms), fin(bms), hv.n > 500 ? "true" : "false");
    for (size_t i = 0; i < hv.n && i < 500; i++) {
        const Hit *h = &hv.a[i];
        char site[SITE_LEN + 1], mask[GUIDE_LEN + 1];
        if (h->strand == '+') memcpy(site, G.seq + h->pos, SITE_LEN);
        else revcomp(G.seq + h->pos, SITE_LEN, site);
        site[SITE_LEN] = 0;
        for (int j = 0; j < GUIDE_LEN; j++) mask[j] = ((h->mask >> j) & 1) ? '1' : '0';
        mask[GUIDE_LEN] = 0;
        bprintf(o, "%s{\"pos\":%llu,\"strand\":\"%c\",\"mm\":%d,\"site\":\"%s\",\"mask\":\"%s\"}", i ? "," : "", LLU(h->pos),
                h->strand, h->mm, site, mask);
    }
    bputs(o, "]}");
    hitvec_free(&hv);
    return 0;
}

static int api_offtarget(const char *q, Buf *o, int *st) {
    if (need_genome(o, st)) return 0;
    long long k, n;
    if (qint(q, "k", 3, 0, 5, &k)) return bad_param(o, st, "k", 0, 5);
    if (qint(q, "n", 10, 1, 50, &n)) return bad_param(o, st, "n", 1, 50);
    TargetVec tv = {0};
    find_targets(&G, kmp_search, &tv);
    if (tv.n == 0) { targetvec_free(&tv); return jerr(o, st, "This genome has no NGG sites to test."); }
    uint64_t *codes = xmalloc((size_t)n * sizeof *codes);
    int nc = 0;
    if (has_si) guide_encode(SI.guide, &codes[nc++]);
    Rng r;
    rng_seed(&r, 1234);
    while (nc < n) {
        Target t = tv.a[rng_below(&r, (uint32_t)tv.n)];
        codes[nc++] = site_code(G.seq, t.pos, t.strand);
    }
    targetvec_free(&tv);
    if (ix_k >= 0) { ot_free(&IX); ix_k = -1; } /* rebuild to time the build honestly */
    double bms;
    const OtIndex *ix = get_index((int)k, &bms);
    Buf rows = {0};
    double tb = 0, ti = 0;
    int same_count = 0;
    for (int i = 0; i < nc; i++) {
        HitVec a = {0}, b = {0};
        double t0 = now_ms();
        ot_brute(&G, codes[i], (int)k, &a);
        double t1 = now_ms();
        ot_query(ix, codes[i], &b);
        double t2 = now_ms();
        hitvec_sort(&a);
        hitvec_sort(&b);
        int same = a.n == b.n && hitvec_signature(&a) == hitvec_signature(&b);
        same_count += same;
        tb += t1 - t0;
        ti += t2 - t1;
        char gs[GUIDE_LEN + 1];
        guide_decode(codes[i], gs);
        bprintf(&rows, "%s{\"guide\":\"%s\",\"hits\":%llu,\"brute_ms\":%.3f,\"index_ms\":%.4f,\"same\":%s}", i ? "," : "", gs,
                LLU(b.n), fin(t1 - t0), fin(t2 - t1), same ? "true" : "false");
        hitvec_free(&a);
        hitvec_free(&b);
    }
    bprintf(o, "{\"k\":%lld,\"n\":%d,\"sites\":%llu,\"build_ms\":%.2f,\"brute_total_ms\":%.2f,\"index_total_ms\":%.3f,"
               "\"brute_per_query_ms\":%.3f,\"index_per_query_ms\":%.4f,\"speedup_query\":%.1f,\"speedup_with_build\":%.2f,"
               "\"identical\":%d,\"of\":%d,\"guides\":[%s]}",
            k, nc, LLU(ix->nsites), fin(bms), fin(tb), fin(ti), fin(tb / nc), fin(ti / nc), fin(ratio(tb, ti)),
            fin(ratio(tb, bms + ti)), same_count, nc, rows.s ? rows.s : "");
    free(rows.s);
    free(codes);
    return 0;
}

static int api_cpg(Buf *o, int *st) {
    if (need_genome(o, st)) return 0;
    IslandVec iv = {0};
    double t0 = now_ms();
    cpg_islands(G.seq, G.n, &iv);
    double ms = now_ms() - t0;
    size_t bins = G.n < 200 ? 1 : 200;
    double *gc = xmalloc(bins * sizeof *gc);
    gc_track(G.seq, G.n, bins, gc);
    bprintf(o, "{\"n\":%llu,\"count\":%llu,\"ms\":%.3f,\"bins\":[", LLU(G.n), LLU(iv.n), fin(ms));
    for (size_t i = 0; i < bins; i++) bprintf(o, "%s%.3f", i ? "," : "", gc[i]);
    bputs(o, "],\"islands\":[");
    for (size_t i = 0; i < iv.n && i < 500; i++)
        bprintf(o, "%s{\"start\":%llu,\"end\":%llu,\"gc\":%.3f,\"oe\":%.2f}", i ? "," : "", LLU(iv.a[i].start), LLU(iv.a[i].end),
                fin(iv.a[i].gc), fin(iv.a[i].oe));
    bputs(o, "]}");
    free(gc);
    islandvec_free(&iv);
    return 0;
}

/* b = a with about 6% substitutions, 3% deletions, 3% insertions, cut to at most la letters */
static char *make_variant(const char *a, size_t la, Rng *r) {
    static const char L[4] = {'A', 'C', 'G', 'T'};
    char *b = xmalloc(la * 2 + 2);
    size_t n = 0;
    for (size_t i = 0; i < la && n < la; i++) {
        double u = rng_unit(r);
        if (u < 0.03) continue;
        if (u < 0.06) b[n++] = L[rng_below(r, 4)];
        b[n++] = u < 0.12 && u >= 0.06 ? L[(base2(a[i]) + 1 + (int)rng_below(r, 3)) % 4] : a[i];
    }
    if (n > la) n = la;
    if (n == 0) b[n++] = 'A';
    b[n] = 0;
    return b;
}

static int clean_seq(const char *in, char *out, size_t cap) { /* returns length or -1 */
    size_t n = 0;
    for (; *in; in++) {
        if (*in == ' ' || *in == '\n' || *in == '\r') continue;
        if (base2(*in) < 0 || n + 1 >= cap) return -1;
        out[n++] = (char)(*in & ~32);
    }
    out[n] = 0;
    return n ? (int)n : -1;
}

static int api_align(const char *q, Buf *o, int *st) {
    char *a, *b;
    size_t m, n;
    char qa[8000], qb[8000];
    int has_a = qstr(q, "a", qa, sizeof qa), has_b = qstr(q, "b", qb, sizeof qb);
    if (has_a != has_b) return jerr(o, st, "Give both sequences: a and b.");
    if (has_a) {
        char ca[3001], cb[3001];
        int la = clean_seq(qa, ca, sizeof ca), lb = clean_seq(qb, cb, sizeof cb);
        if (la < 0 || lb < 0)
            return jerr(o, st, "Each sequence must have 1 to 3000 letters, all A, C, G or T.");
        m = (size_t)la; n = (size_t)lb;
        a = xmalloc(m + 1); b = xmalloc(n + 1);
        memcpy(a, ca, m + 1); memcpy(b, cb, n + 1);
    } else {
        long long len;
        if (qint(q, "len", 1000, 50, 8000, &len)) return bad_param(o, st, "len", 50, 8000);
        m = (size_t)len;
        a = xmalloc(m + 1);
        Rng r;
        rng_seed(&r, 31337);
        if (loaded && G.n >= m + 10) memcpy(a, G.seq + G.n / 2 - (m / 2 < G.n / 2 ? m / 2 : 0), m);
        else for (size_t i = 0; i < m; i++) a[i] = "ACGT"[rng_below(&r, 4)];
        a[m] = 0;
        b = make_variant(a, m, &r);
        n = strlen(b);
    }
    Alignment nw = {0}, hb;
    double t0 = now_ms();
    int nw_ok = nw_align(a, m, b, n, &nw) == 0;
    double t_nw = now_ms() - t0;
    t0 = now_ms();
    hirschberg_align(a, m, b, n, &hb);
    double t_hb = now_ms() - t0;
    size_t swb;
    t0 = now_ms();
    int sw = sw_score(a, m, b, n, &swb);
    double t_sw = now_ms() - t0;
    size_t mt = 0, mis = 0, gaps = 0;
    for (size_t i = 0; i < hb.len; i++) {
        if (hb.a[i] == '-' || hb.b[i] == '-') gaps++;
        else if (hb.a[i] == hb.b[i]) mt++;
        else mis++;
    }
    size_t shown = hb.len < 600 ? hb.len : 600;
    bprintf(o, "{\"m\":%llu,\"n\":%llu,\"hirschberg\":{\"score\":%d,\"bytes\":%llu,\"ms\":%.2f},\"sw\":{\"score\":%d,\"bytes\":%llu,\"ms\":%.2f},",
            LLU(m), LLU(n), hb.score, LLU(hb.bytes), fin(t_hb), sw, LLU(swb), fin(t_sw));
    if (nw_ok) bprintf(o, "\"nw\":{\"score\":%d,\"bytes\":%llu,\"ms\":%.2f},\"equal\":%s,", nw.score, LLU(nw.bytes), fin(t_nw), nw.score == hb.score ? "true" : "false");
    else bputs(o, "\"nw\":null,\"equal\":null,");
    bprintf(o, "\"matches\":%llu,\"mismatches\":%llu,\"gaps\":%llu,\"length\":%llu,\"shown\":%llu,\"a\":\"%.*s\",\"b\":\"%.*s\"}",
            LLU(mt), LLU(mis), LLU(gaps), LLU(hb.len), LLU(shown), (int)shown, hb.a, (int)shown, hb.b);
    if (nw_ok) alignment_free(&nw);
    alignment_free(&hb);
    free(a);
    free(b);
    return 0;
}

static int api_vault(const char *q, Buf *o, int *st) {
    if (need_genome(o, st)) return 0;
    long long bits;
    if (qint(q, "bits", 64, 40, 64, &bits)) return bad_param(o, st, "bits", 40, 64);
    const char *f0 = "crispr_tmp_plain.fa", *f1 = "crispr_tmp_enc.bin", *f2 = "crispr_tmp_back.fa", *f3 = "crispr_tmp_bad.fa";
    char err[256], e2[256], e3[256];
    Genome part = {0};
    part.seq = G.seq;
    part.n = G.n < 200000 ? G.n : 200000;
    snprintf(part.name, sizeof part.name, "vault_sample");
    Rng r;
    rng_seed(&r, os_random64());
    RsaKey k1, k2;
    rsa_keygen((int)bits, &r, &k1);
    rsa_keygen((int)bits, &r, &k2);
    VaultStats stt = {0};
    int ok = genome_save(f0, &part, err, sizeof err) == 0 && vault_encrypt_file(f0, f1, &k1, &r, &stt, err, sizeof err) == 0;
    int same = 0, rej_other = 0, rej_exp = 0;
    e2[0] = e3[0] = 0;
    if (ok) {
        same = vault_decrypt_file(f1, f2, &k1, err, sizeof err) == 0 && files_equal(f0, f2);
        rej_other = vault_decrypt_file(f1, f3, &k2, e2, sizeof e2) != 0;
        RsaKey bad = k1;
        bad.d ^= 0x10;
        rej_exp = vault_decrypt_file(f1, f3, &bad, e3, sizeof e3) != 0;
    }
    remove(f0); remove(f1); remove(f2); remove(f3);
    if (!ok) return jerr(o, st, err);
    bprintf(o, "{\"bits\":%lld,\"modulus\":\"%llu\",\"e\":%llu,\"original\":%llu,\"compressed\":%llu,\"encrypted\":%llu,\"restored\":%llu,"
               "\"checks\":[{\"name\":\"Restored file is byte-identical to the original\",\"pass\":%s},"
               "{\"name\":\"A different key pair is rejected\",\"pass\":%s,\"message\":",
            bits, LLU(k1.n), LLU(k1.e), LLU(stt.orig), LLU(stt.comp), LLU(stt.enc), LLU(same ? stt.orig : 0), same ? "true" : "false", rej_other ? "true" : "false");
    bstr(o, e2);
    bprintf(o, "},{\"name\":\"Right modulus, wrong private exponent is rejected\",\"pass\":%s,\"message\":", rej_exp ? "true" : "false");
    bstr(o, e3);
    bputs(o, "}]}");
    return 0;
}

static int api_attack(Buf *o, int *st) {
    (void)st;
    AttackRow rows[16];
    int n = attack_table(12, 32, 2, 30, 2024, rows, 16);
    bputs(o, "{\"rows\":[");
    for (int i = 0; i < n; i++)
        bprintf(o, "%s{\"prime_bits\":%d,\"modulus_bits\":%d,\"avg_steps\":%.1f,\"theory\":%.1f,\"avg_ms\":%.4f,\"recovered\":%d,\"trials\":%d}",
                i ? "," : "", rows[i].prime_bits, rows[i].prime_bits * 2, fin(rows[i].avg_steps), fin(rows[i].theory), fin(rows[i].avg_ms), rows[i].ok, rows[i].trials);
    bputs(o, "],\"extrapolation\":[");
    const AttackRow *last = &rows[n - 1];
    double rate = ratio(last->avg_steps, last->avg_ms) * 1000.0; /* steps per second, measured */
    int targets[4] = {64, 128, 512, 1024};
    for (int i = 0; i < 4; i++) {
        double log2s = log2(last->avg_steps) + (targets[i] - last->prime_bits) / 2.0;
        double log10_years = log2s * log10(2.0) - log10(rate) - log10(31557600.0);
        bprintf(o, "%s{\"prime_bits\":%d,\"modulus_bits\":%d,\"log2_steps\":%.1f,\"log10_years\":%.1f}", i ? "," : "", targets[i], targets[i] * 2, log2s, log10_years);
    }
    bprintf(o, "],\"steps_per_second\":%.0f}", fin(rate));
    return 0;
}

static int api_bench(Buf *o, int *st) {
    size_t sizes[4] = {100000, 300000, 1000000, 4600000};
    size_t lens[4] = {500, 1000, 2000, 4000};
    GenomeBenchRow gr[4];
    AlignBenchRow ar[4];
    if (bench_genomes(sizes, 4, gr) || bench_align(lens, 4, ar)) return jerr(o, st, "The benchmark failed.");
    bputs(o, "{\"genomes\":[");
    for (int i = 0; i < 4; i++)
        bprintf(o, "%s{\"n\":%llu,\"sites\":%llu,\"kmp_ms\":%.3f,\"bm_ms\":%.3f,\"rk_ms\":%.3f,\"build_ms\":%.2f,\"query_ms\":%.4f,\"brute_ms\":%.3f,\"cpg_ms\":%.3f}",
                i ? "," : "", LLU(gr[i].n), LLU(gr[i].sites), fin(gr[i].kmp_ms), fin(gr[i].bm_ms), fin(gr[i].rk_ms), fin(gr[i].build_ms), fin(gr[i].query_ms), fin(gr[i].brute_ms), fin(gr[i].cpg_ms));
    bputs(o, "],\"align\":[");
    for (int i = 0; i < 4; i++)
        bprintf(o, "%s{\"len\":%llu,\"nw_ms\":%.2f,\"hb_ms\":%.2f,\"nw_bytes\":%llu,\"hb_bytes\":%llu}", i ? "," : "", LLU(ar[i].len), fin(ar[i].nw_ms), fin(ar[i].hb_ms), LLU(ar[i].nw_bytes), LLU(ar[i].hb_bytes));
    bputs(o, "]}");
    return 0;
}

/* ---------- HTTP plumbing ---------- */
static const char *STATUS_TEXT(int s) {
    switch (s) {
    case 200: return "OK";
    case 400: return "Bad Request";
    case 404: return "Not Found";
    case 413: return "Payload Too Large";
    case 431: return "Request Header Fields Too Large";
    default: return "Internal Server Error";
    }
}

static void send_all(sock_t s, const char *p, size_t n) {
    while (n) {
        int w = (int)send(s, p, n > 65536 ? 65536 : (int)n, 0);
        if (w <= 0) return;
        p += w;
        n -= (size_t)w;
    }
}

static void respond(sock_t s, int status, const char *ctype, const char *body, size_t n) {
    char h[256];
    int l = snprintf(h, sizeof h,
                     "HTTP/1.1 %d %s\r\nContent-Type: %s\r\nContent-Length: %llu\r\nConnection: close\r\nCache-Control: no-store\r\n\r\n",
                     status, STATUS_TEXT(status), ctype, LLU(n));
    send_all(s, h, (size_t)l);
    send_all(s, body, n);
}

static void respond_json(sock_t s, int status, const char *msg) {
    Buf b = {0};
    int st = status;
    jerr(&b, &st, msg);
    respond(s, status, "application/json", b.s, b.n);
    free(b.s);
}

static int ci_prefix(const char *s, const char *p) { /* p is lowercase */
    for (; *p; s++, p++) {
        char a = *s;
        if (a >= 'A' && a <= 'Z') a = (char)(a + 32);
        if (a != *p) return 0;
    }
    return 1;
}

static long content_length(const char *hdr, size_t hlen) {
    for (size_t i = 0; i + 15 < hlen; i++)
        if ((i == 0 || hdr[i - 1] == '\n') && ci_prefix(hdr + i, "content-length:")) return atol(hdr + i + 15);
    return 0;
}

static char g_webdir[512] = "web";

static void serve_index(sock_t s) {
    char path[600];
    snprintf(path, sizeof path, "%s/index.html", g_webdir);
    FILE *f = fopen(path, "rb");
    if (!f) {
        const char *m = "web/index.html was not found. Start the server from the crispr_vault folder.";
        respond(s, 500, "text/plain; charset=utf-8", m, strlen(m));
        return;
    }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buf = xmalloc((size_t)sz + 1);
    size_t got = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    respond(s, 200, "text/html; charset=utf-8", buf, got);
    free(buf);
}

static void handle(sock_t s, const char *method, const char *path, const char *query, const char *body, size_t blen) {
    if (!strcmp(path, "/") || !strcmp(path, "/index.html")) { serve_index(s); return; }
    Buf o = {0};
    int st = 200;
    if (!strcmp(path, "/api/status")) status_json(&o);
    else if (!strcmp(path, "/api/synth")) api_synth(query, &o, &st);
    else if (!strcmp(path, "/api/upload")) {
        if (strcmp(method, "POST")) jerr(&o, &st, "Upload needs a POST request with the FASTA file as the body.");
        else api_upload(query, body, blen, &o, &st);
    }
    else if (!strcmp(path, "/api/find")) api_find(&o, &st);
    else if (!strcmp(path, "/api/rank")) api_rank(query, &o, &st);
    else if (!strcmp(path, "/api/hits")) api_hits(query, &o, &st);
    else if (!strcmp(path, "/api/offtarget")) api_offtarget(query, &o, &st);
    else if (!strcmp(path, "/api/cpg")) api_cpg(&o, &st);
    else if (!strcmp(path, "/api/align")) api_align(query, &o, &st);
    else if (!strcmp(path, "/api/vault")) api_vault(query, &o, &st);
    else if (!strcmp(path, "/api/attack")) api_attack(&o, &st);
    else if (!strcmp(path, "/api/bench")) api_bench(&o, &st);
    else { respond_json(s, 404, "No such page or endpoint."); return; }
    respond(s, st, "application/json", o.s ? o.s : "{}", o.n);
    free(o.s);
}

static void serve_one(sock_t c) {
    size_t cap = 16384, len = 0;
    char *buf = xmalloc(cap);
    char *hend = NULL;
    for (;;) {
        if (cap - len < 2) { cap *= 2; buf = realloc(buf, cap); if (!buf) exit(1); }
        int r = (int)recv(c, buf + len, (int)(cap - len - 1), 0);
        if (r <= 0) break;
        len += (size_t)r;
        buf[len] = 0;
        hend = strstr(buf, "\r\n\r\n");
        if (hend) break;
        if (len >= 65536) { respond_json(c, 431, "The request headers are too large."); free(buf); return; }
    }
    if (!hend) { free(buf); return; }
    size_t hlen = (size_t)(hend - buf) + 4;
    long cl = content_length(buf, hlen);
    if (cl < 0 || (unsigned long)cl > MAX_BODY) { respond_json(c, 413, "The file is larger than the 64 MB limit."); free(buf); return; }
    size_t total = hlen + (size_t)cl;
    if (cap < total + 1) { cap = total + 1; buf = realloc(buf, cap); if (!buf) exit(1); }
    while (len < total) {
        int r = (int)recv(c, buf + len, (int)(total - len > 1048576 ? 1048576 : total - len), 0);
        if (r <= 0) break;
        len += (size_t)r;
    }
    if (len < total) { respond_json(c, 400, "The upload ended before the whole file arrived."); free(buf); return; }
    buf[len] = 0;
    char *method = buf, *sp = strchr(buf, ' ');
    if (!sp) { respond_json(c, 400, "Malformed request."); free(buf); return; }
    *sp = 0;
    char *path = sp + 1, *sp2 = strchr(path, ' ');
    if (!sp2) { respond_json(c, 400, "Malformed request."); free(buf); return; }
    *sp2 = 0;
    const char *query = "";
    char *qm = strchr(path, '?');
    if (qm) { *qm = 0; query = qm + 1; }
    handle(c, method, path, query, buf + hlen, (size_t)cl);
    free(buf);
}

int serve_main(int port, const char *webdir) {
    snprintf(g_webdir, sizeof g_webdir, "%s", webdir);
#ifdef _WIN32
    WSADATA w;
    if (WSAStartup(MAKEWORD(2, 2), &w) != 0) { fprintf(stderr, "error: Winsock could not start.\n"); return 1; }
#else
    signal(SIGPIPE, SIG_IGN);
#endif
    sock_t ls = socket(AF_INET, SOCK_STREAM, 0);
    if (ls == BADSOCK) { fprintf(stderr, "error: cannot create a socket.\n"); return 1; }
#ifndef _WIN32
    int yes = 1;
    setsockopt(ls, SOL_SOCKET, SO_REUSEADDR, (const char *)&yes, sizeof yes);
#endif
    struct sockaddr_in a;
    memset(&a, 0, sizeof a);
    a.sin_family = AF_INET;
    a.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    a.sin_port = htons((unsigned short)port);
    if (bind(ls, (struct sockaddr *)&a, sizeof a) != 0 || listen(ls, 8) != 0) {
        fprintf(stderr, "error: cannot listen on 127.0.0.1:%d (is the port already in use?).\n", port);
        return 1;
    }
    printf("Open http://127.0.0.1:%d/ in a browser. Press Ctrl+C to stop.\n", port);
    fflush(stdout);
    for (;;) {
        sock_t c = accept(ls, NULL, NULL);
        if (c == BADSOCK) continue;
        serve_one(c);
        CLOSESOCK(c);
    }
}
