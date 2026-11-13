/* vault.c: compress, encrypt, decrypt, and break (a toy version of) the lock.
 *
 * TEACHING IMPLEMENTATION. 64-bit RSA modulus, a simple xorshift keystream,
 * not for real data.
 *
 * Huffman coding   O(n + sigma log sigma). The tree is rebuilt from a stored
 *                  frequency table with fixed tie-breaking, so the decoder
 *                  builds exactly the tree the encoder used.
 * Miller-Rabin     O(log^2 n) per base with 12 fixed bases (exact below 2^64).
 * powmod           O(log e) multiplications (square and multiply, __int128).
 * Extended Euclid  O(log n) steps, gives the private exponent d = e^-1 mod phi.
 * File layout      "GVLT" | modulus | RSA(session seed) | payload length | payload
 *                  payload = keystream XOR (huffman blob | FNV-1a of original)
 * Pollard's rho    about O(sqrt(p)) steps for the smallest prime factor p,
 *                  i.e. 2^(bits/2) for a bits-bit prime. Breaks the key by
 *                  factoring the modulus.
 */
#ifdef _WIN32
#define _CRT_RAND_S
#endif
#include "crispr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------- number theory ---------- */
uint64_t os_random64(void) {
#ifdef _WIN32
    unsigned int a = 0, b = 0;
    rand_s(&a);
    rand_s(&b);
    return ((uint64_t)a << 32) | b;
#else
    uint64_t v = 0;
    FILE *f = fopen("/dev/urandom", "rb");
    if (f) { if (fread(&v, sizeof v, 1, f) != 1) v = 0; fclose(f); }
    if (!v) v = (uint64_t)(now_ms() * 1000.0) * 2654435761ULL;
    return v;
#endif
}

uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) {
    return (uint64_t)((unsigned __int128)a * b % m);
}

uint64_t powmod(uint64_t b, uint64_t e, uint64_t m) {
    uint64_t r = 1 % m;
    b %= m;
    while (e) {
        if (e & 1) r = mulmod(r, b, m);
        b = mulmod(b, b, m);
        e >>= 1;
    }
    return r;
}

int is_prime64(uint64_t n) {
    static const uint64_t B[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return 0;
    for (int i = 0; i < 12; i++) {
        if (n == B[i]) return 1;
        if (n % B[i] == 0) return 0;
    }
    uint64_t d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        uint64_t x = powmod(B[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int comp = 1;
        for (int r = 1; r < s; r++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { comp = 0; break; }
        }
        if (comp) return 0;
    }
    return 1;
}

static uint64_t gcd64(uint64_t a, uint64_t b) { while (b) { uint64_t t = a % b; a = b; b = t; } return a; }

static int modinv(uint64_t a, uint64_t m, uint64_t *out) { /* extended Euclid */
    __int128 t = 0, nt = 1, r = (__int128)m, nr = (__int128)a;
    while (nr != 0) {
        __int128 q = r / nr, tmp = t - q * nt;
        t = nt; nt = tmp;
        tmp = r - q * nr;
        r = nr; nr = tmp;
    }
    if (r != 1) return -1;
    if (t < 0) t += (__int128)m;
    *out = (uint64_t)t;
    return 0;
}

static uint64_t random_prime(int bits, Rng *r) {
    uint64_t top = 1ULL << (bits - 1), mask = top - 1;
    for (;;) {
        uint64_t x = top | (rng_next(r) & mask) | 1;
        if (is_prime64(x)) return x;
    }
}

int rsa_keygen(int modbits, Rng *r, RsaKey *k) {
    if (modbits < 16 || modbits > 64) return -1;
    int pb = modbits / 2;
    memset(k, 0, sizeof *k);
    for (;;) {
        uint64_t p = random_prime(pb, r), q = random_prime(pb, r);
        if (p == q) continue;
        uint64_t phi = (p - 1) * (q - 1);
        static const uint64_t ES[5] = {65537, 257, 17, 5, 3};
        uint64_t e = 0, d = 0;
        for (int i = 0; i < 5; i++)
            if (ES[i] < phi && gcd64(ES[i], phi) == 1) { e = ES[i]; break; }
        if (!e || modinv(e, phi, &d) != 0) continue;
        k->p = p; k->q = q; k->n = p * q; k->e = e; k->d = d; k->bits = modbits;
        return 0;
    }
}

int rsa_key_save(const char *path, const RsaKey *k, int priv, char *err, size_t errsz) {
    FILE *f = fopen(path, "wb");
    if (!f) { if (err) snprintf(err, errsz, "Cannot write key file '%s'.", path); return -1; }
    fprintf(f, "%s %llu %llu\n", priv ? "RSA-PRIV" : "RSA-PUB", (unsigned long long)k->n,
            (unsigned long long)(priv ? k->d : k->e));
    fclose(f);
    return 0;
}

int rsa_key_load(const char *path, RsaKey *k, char *err, size_t errsz) {
    FILE *f = fopen(path, "rb");
    char tag[16];
    unsigned long long n, x;
    memset(k, 0, sizeof *k);
    if (!f) { if (err) snprintf(err, errsz, "Cannot open key file '%s'.", path); return -1; }
    int got = fscanf(f, "%15s %llu %llu", tag, &n, &x);
    fclose(f);
    if (got != 3 || (strcmp(tag, "RSA-PUB") && strcmp(tag, "RSA-PRIV")) || n < 2) {
        if (err) snprintf(err, errsz, "'%s' is not a key file made by this program.", path);
        return -1;
    }
    k->n = n;
    if (!strcmp(tag, "RSA-PUB")) k->e = x; else k->d = x;
    return 0;
}

uint64_t fnv1a(const uint8_t *p, size_t n) {
    uint64_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}

/* ---------- Huffman ---------- */
typedef struct { uint64_t freq; int left, right, parent; uint8_t sym; } HNode;

static int hcmp(const HNode *nd, int a, int b) { /* heap order: freq, then index */
    if (nd[a].freq != nd[b].freq) return nd[a].freq < nd[b].freq ? -1 : 1;
    return a < b ? -1 : (a > b);
}
static void hpush(int *h, int *sz, int x, const HNode *nd) {
    int c = (*sz)++;
    h[c] = x;
    while (c > 0) {
        int p = (c - 1) / 2;
        if (hcmp(nd, h[c], h[p]) < 0) { int t = h[c]; h[c] = h[p]; h[p] = t; c = p; } else break;
    }
}
static int hpop(int *h, int *sz, const HNode *nd) {
    int top = h[0];
    h[0] = h[--(*sz)];
    int i = 0;
    for (;;) {
        int l = 2 * i + 1, r = l + 1, s = i;
        if (l < *sz && hcmp(nd, h[l], h[s]) < 0) s = l;
        if (r < *sz && hcmp(nd, h[r], h[s]) < 0) s = r;
        if (s == i) break;
        int t = h[i]; h[i] = h[s]; h[s] = t; i = s;
    }
    return top;
}

/* leaves 0..nsym-1 in symbol order, internal nodes after; returns root */
static int build_tree(const uint32_t *freq, HNode *nd, int *nsym_out) {
    int nsym = 0;
    for (int s = 0; s < 256; s++)
        if (freq[s]) { nd[nsym] = (HNode){freq[s], -1, -1, -1, (uint8_t)s}; nsym++; }
    *nsym_out = nsym;
    if (nsym == 0) return -1;
    int heap[512], hs = 0, total = nsym;
    for (int i = 0; i < nsym; i++) hpush(heap, &hs, i, nd);
    while (hs > 1) {
        int a = hpop(heap, &hs, nd), b = hpop(heap, &hs, nd);
        nd[total] = (HNode){nd[a].freq + nd[b].freq, a, b, -1, 0};
        nd[a].parent = nd[b].parent = total;
        hpush(heap, &hs, total, nd);
        total++;
    }
    return heap[0];
}

static void put_u64(uint8_t *p, uint64_t v) { for (int i = 0; i < 8; i++) p[i] = (uint8_t)(v >> (8 * i)); }
static uint64_t get_u64(const uint8_t *p) { uint64_t v = 0; for (int i = 0; i < 8; i++) v |= (uint64_t)p[i] << (8 * i); return v; }

int huff_compress(const uint8_t *in, size_t n, uint8_t **out, size_t *outn) {
    if ((uint64_t)n > 0xFFFFFFFFULL) return -1;
    uint32_t freq[256] = {0};
    for (size_t i = 0; i < n; i++) freq[in[i]]++;
    HNode nd[512];
    int nsym, root = build_tree(freq, nd, &nsym);
    static uint8_t bits[256][256];
    int len[256] = {0}, leaf_of[256];
    if (nsym > 1) {
        for (int i = 0; i < nsym; i++) {
            uint8_t tmp[256];
            int l = 0, x = i;
            while (nd[x].parent >= 0) {
                int par = nd[x].parent;
                tmp[l++] = (nd[par].right == x);
                x = par;
            }
            int s = nd[i].sym;
            len[s] = l;
            for (int j = 0; j < l; j++) bits[s][j] = tmp[l - 1 - j];
            leaf_of[s] = i;
        }
    }
    (void)root;
    uint64_t nbits = 0;
    if (nsym > 1) for (int s = 0; s < 256; s++) nbits += (uint64_t)freq[s] * (uint64_t)len[s];
    size_t hdr = 3 + 8 + 2 + (size_t)nsym * 5;
    size_t total = hdr + (size_t)((nbits + 7) / 8);
    uint8_t *o = xcalloc(total, 1);
    memcpy(o, "HF1", 3);
    put_u64(o + 3, n);
    o[11] = (uint8_t)(nsym & 255); o[12] = (uint8_t)(nsym >> 8);
    size_t p = 13;
    for (int s = 0; s < 256; s++)
        if (freq[s]) {
            o[p++] = (uint8_t)s;
            for (int j = 0; j < 4; j++) o[p++] = (uint8_t)(freq[s] >> (8 * j));
        }
    if (nsym > 1) {
        uint64_t bp = 0;
        for (size_t i = 0; i < n; i++) {
            int s = in[i];
            for (int j = 0; j < len[s]; j++, bp++)
                if (bits[s][j]) o[hdr + (size_t)(bp >> 3)] |= (uint8_t)(0x80 >> (bp & 7));
        }
    }
    (void)leaf_of;
    *out = o;
    *outn = total;
    return 0;
}

int huff_decompress(const uint8_t *in, size_t n, uint8_t **out, size_t *outn) {
    if (n < 13 || memcmp(in, "HF1", 3) != 0) return -1;
    uint64_t cnt = get_u64(in + 3);
    int nsym = in[11] | (in[12] << 8);
    if (nsym > 256 || (cnt > 0 && nsym == 0) || n < 13 + (size_t)nsym * 5) return -1;
    uint32_t freq[256] = {0};
    size_t p = 13;
    uint64_t sum = 0;
    int prev = -1;
    for (int i = 0; i < nsym; i++) {
        int s = in[p++];
        if (s <= prev) return -1;
        prev = s;
        uint32_t f = 0;
        for (int j = 0; j < 4; j++) f |= (uint32_t)in[p++] << (8 * j);
        if (f == 0) return -1;
        freq[s] = f;
        sum += f;
    }
    if (sum != cnt || cnt > (1ULL << 30)) return -1;
    size_t hdr = p;
    uint64_t avail = (uint64_t)(n - hdr) * 8;
    if (nsym > 1 && cnt > avail) return -1;
    uint8_t *o = xmalloc((size_t)cnt + 1);
    if (cnt == 0) { *out = o; *outn = 0; return 0; }
    HNode nd[512];
    int ns2, root = build_tree(freq, nd, &ns2);
    if (nsym == 1) {
        memset(o, nd[0].sym, (size_t)cnt);
    } else {
        uint64_t bp = 0;
        for (uint64_t i = 0; i < cnt; i++) {
            int x = root;
            while (nd[x].left >= 0) {
                if (bp >= avail) { free(o); return -1; }
                int bit = (in[hdr + (size_t)(bp >> 3)] >> (7 - (bp & 7))) & 1;
                bp++;
                x = bit ? nd[x].right : nd[x].left;
            }
            o[i] = nd[x].sym;
        }
    }
    *out = o;
    *outn = (size_t)cnt;
    return 0;
}

/* ---------- stream cipher and container ---------- */
static void keystream_xor(uint8_t *buf, size_t n, uint32_t seed) {
    uint64_t x = mix64((uint64_t)seed ^ 0xA5A5A5A55A5A5A5AULL);
    if (!x) x = 1;
    for (size_t i = 0; i < n; i++) {
        x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
        buf[i] ^= (uint8_t)((x * 0x2545F4914F6CDD1DULL) >> 56);
    }
}

int vault_encrypt(const uint8_t *in, size_t n, const RsaKey *k, Rng *r,
                  uint8_t **out, size_t *outn, VaultStats *st) {
    if (k->n <= 0xFFFFFFFFULL || !k->e) return -1; /* seed must fit below the modulus */
    uint8_t *comp;
    size_t cn;
    if (huff_compress(in, n, &comp, &cn) != 0) return -1;
    size_t plen = cn + 8;
    uint8_t *o = xmalloc(28 + plen);
    memcpy(o, "GVLT", 4);
    put_u64(o + 4, k->n);
    uint32_t seed = (uint32_t)rng_next(r);
    put_u64(o + 12, powmod(seed, k->e, k->n));
    put_u64(o + 20, plen);
    memcpy(o + 28, comp, cn);
    put_u64(o + 28 + cn, fnv1a(in, n));
    keystream_xor(o + 28, plen, seed);
    free(comp);
    *out = o;
    *outn = 28 + plen;
    if (st) { st->orig = n; st->comp = cn; st->enc = 28 + plen; }
    return 0;
}

static void say(char *err, size_t errsz, const char *m) { if (err && errsz) snprintf(err, errsz, "%s", m); }

int vault_decrypt(const uint8_t *in, size_t n, const RsaKey *k,
                  uint8_t **out, size_t *outn, char *err, size_t errsz) {
    if (n < 36 || memcmp(in, "GVLT", 4) != 0) { say(err, errsz, "This is not a vault file."); return -1; }
    if (get_u64(in + 4) != k->n) { say(err, errsz, "Wrong key: the modulus does not match this file."); return -1; }
    uint64_t plen = get_u64(in + 20);
    if (plen != n - 28) { say(err, errsz, "The vault file is damaged (length mismatch)."); return -1; }
    uint64_t seed = powmod(get_u64(in + 12), k->d, k->n);
    if (seed > 0xFFFFFFFFULL) { say(err, errsz, "Wrong key: the session seed did not decrypt."); return -1; }
    uint8_t *pl = xmalloc((size_t)plen);
    memcpy(pl, in + 28, (size_t)plen);
    keystream_xor(pl, (size_t)plen, (uint32_t)seed);
    uint8_t *data;
    size_t dn;
    if (huff_decompress(pl, (size_t)plen - 8, &data, &dn) != 0) {
        free(pl);
        say(err, errsz, "Wrong key or damaged file: the decrypted data is not valid.");
        return -1;
    }
    uint64_t want = get_u64(pl + plen - 8);
    free(pl);
    if (fnv1a(data, dn) != want) {
        free(data);
        say(err, errsz, "Wrong key or damaged file: the checksum does not match.");
        return -1;
    }
    *out = data;
    *outn = dn;
    return 0;
}

/* ---------- files ---------- */
static int read_all(const char *path, uint8_t **buf, size_t *n, char *err, size_t errsz) {
    FILE *f = fopen(path, "rb");
    if (!f) { if (err) snprintf(err, errsz, "Cannot open file '%s'.", path); return -1; }
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    *buf = xmalloc(sz > 0 ? (size_t)sz : 1);
    *n = fread(*buf, 1, sz > 0 ? (size_t)sz : 0, f);
    fclose(f);
    return 0;
}
static int write_all(const char *path, const uint8_t *buf, size_t n, char *err, size_t errsz) {
    FILE *f = fopen(path, "wb");
    if (!f) { if (err) snprintf(err, errsz, "Cannot write file '%s'.", path); return -1; }
    size_t w = n ? fwrite(buf, 1, n, f) : 0;
    fclose(f);
    if (w != n) { say(err, errsz, "Short write (disk full?)."); return -1; }
    return 0;
}

int vault_encrypt_file(const char *in, const char *out, const RsaKey *k, Rng *r,
                       VaultStats *st, char *err, size_t errsz) {
    uint8_t *b, *o;
    size_t n, on;
    if (read_all(in, &b, &n, err, errsz) != 0) return -1;
    if (vault_encrypt(b, n, k, r, &o, &on, st) != 0) {
        free(b);
        say(err, errsz, "Encryption failed: use a public key with a modulus above 2^32.");
        return -1;
    }
    free(b);
    int rc = write_all(out, o, on, err, errsz);
    free(o);
    return rc;
}

int vault_decrypt_file(const char *in, const char *out, const RsaKey *k, char *err, size_t errsz) {
    uint8_t *b, *o;
    size_t n, on;
    if (read_all(in, &b, &n, err, errsz) != 0) return -1;
    int rc = vault_decrypt(b, n, k, &o, &on, err, errsz);
    free(b);
    if (rc != 0) return -1;
    rc = write_all(out, o, on, err, errsz);
    free(o);
    return rc;
}

int files_equal(const char *p1, const char *p2) {
    FILE *a = fopen(p1, "rb"), *b = fopen(p2, "rb");
    int eq = a && b;
    while (eq) {
        unsigned char x[4096], y[4096];
        size_t na = fread(x, 1, sizeof x, a), nb = fread(y, 1, sizeof y, b);
        if (na != nb || memcmp(x, y, na) != 0) eq = 0;
        if (na == 0) break;
    }
    if (a) fclose(a);
    if (b) fclose(b);
    return eq;
}

/* ---------- Pollard's rho ---------- */
uint64_t pollard_rho(uint64_t n, uint64_t *steps) {
    uint64_t cnt = 0;
    if (n % 2 == 0) { if (steps) *steps = 0; return 2; }
    for (uint64_t c = 1;; c++) {
        uint64_t x = 2, y = 2, d = 1;
        while (d == 1) {
            x = (uint64_t)(((unsigned __int128)mulmod(x, x, n) + c) % n);
            y = (uint64_t)(((unsigned __int128)mulmod(y, y, n) + c) % n);
            y = (uint64_t)(((unsigned __int128)mulmod(y, y, n) + c) % n);
            d = gcd64(x > y ? x - y : y - x, n);
            cnt++;
        }
        if (d != n) { if (steps) *steps = cnt; return d; }
    }
}

int attack_table(int minb, int maxb, int step, int trials, uint64_t seed, AttackRow *rows, int cap) {
    Rng r;
    rng_seed(&r, seed);
    int nr = 0;
    for (int b = minb; b <= maxb && nr < cap; b += step) {
        AttackRow row = {b, trials, 0, 0.0, 0.0, 0.0};
        double steps_sum = 0, ms_sum = 0;
        for (int t = 0; t < trials; t++) {
            RsaKey k;
            rsa_keygen(2 * b, &r, &k);
            uint64_t st;
            double t0 = now_ms();
            uint64_t f = pollard_rho(k.n, &st);
            ms_sum += now_ms() - t0;
            steps_sum += (double)st;
            uint64_t g = k.n / f, d2;
            if (f > 1 && f < k.n && modinv(k.e, (f - 1) * (g - 1), &d2) == 0 && d2 == k.d) row.ok++;
        }
        row.avg_steps = steps_sum / trials;
        row.avg_ms = ms_sum / trials;
        row.theory = (double)(1ULL << (b / 2)) * ((b & 1) ? 1.4142135623730951 : 1.0);
        rows[nr++] = row;
    }
    return nr;
}
