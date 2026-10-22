/* motif_mine.c — short linear motif (6-10 aa) cross-match miner with composition-matched null
 *
 * Goal: test whether EBNA1 contains short human-identical motifs that cannot be produced by
 * shuffling EBNA1 while preserving its local composition (window-preserving shuffle).
 *
 * Method: for each motif length L and mismatch budget d, split the motif into d+1 contiguous
 * pieces (pigeonhole: any pair with <=d mismatches must share at least one piece exactly).
 * Index every human L-window by each piece; look up query windows per piece; verify.
 *
 * Composition gate (optional, --gate 1): a motif must have >=60% distinct residues and no single
 * residue occupying >40% of the motif. This is what separates "order-driven" matches from the
 * low-complexity (G/A-rich, S/P/R-rich) matches that dominate EBNA1.
 *
 * Usage:
 *   motif_mine --query P03211.fasta --db human_swissprot.fasta \
 *              --out-prefix reports/motif --null-b 50 --seed 42 --gate 1 --threads 14
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>

#define MAXSEQ   40000
#define MAXQLEN  200000
#define MAXTHR   64
#define ROWCAP   6000000L   /* cap on detailed real rows per L */

static int8_t CODE[256];
static const char AALET[21] = "ACDEFGHIKLMNPQRSTVWY";

typedef struct { char acc[32]; char gene[32]; int len; uint8_t *s; } Seq;

static Seq   *DB  = NULL;
static int    NDB = 0;
static uint8_t Q[MAXQLEN];
static int    QLEN = 0;
static char   QACC[32], QGENE[32];

static int  GATE = 1;
static int  THREADS = 8;

/* ---------- helpers ---------- */

static void init_code(void) {
    for (int i = 0; i < 256; i++) CODE[i] = -1;
    const char *aa = "ACDEFGHIKLMNPQRSTVWY";
    for (int i = 0; i < 20; i++) { CODE[(unsigned char)aa[i]] = (int8_t)i; CODE[(unsigned char)(aa[i] + 32)] = (int8_t)i; }
}

static char *read_all(const char *path, size_t *n) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
    char *buf = (char *)malloc((size_t)sz + 1);
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) { fprintf(stderr, "read error %s\n", path); exit(1); }
    buf[sz] = 0; fclose(f); *n = (size_t)sz; return buf;
}

static void parse_header(char *h, char *acc, char *gene) {
    acc[0] = 0; gene[0] = 0;
    char *b1 = strchr(h, '|');
    if (b1) {
        char *b2 = strchr(b1 + 1, '|');
        if (b2) { int n = (int)(b2 - b1 - 1); if (n > 31) n = 31; memcpy(acc, b1 + 1, (size_t)n); acc[n] = 0; }
    }
    if (!acc[0]) {
        int n = 0; while (h[n] && h[n] != ' ' && h[n] != '\t' && n < 31) { acc[n] = h[n]; n++; }
        acc[n] = 0;
    }
    char *g = strstr(h, "GN=");
    if (g) {
        g += 3; int n = 0;
        while (g[n] && g[n] != ' ' && g[n] != '\t' && n < 31) { gene[n] = g[n]; n++; }
        gene[n] = 0;
    }
    if (!gene[0]) { int n = 0; while (acc[n]) { gene[n] = acc[n]; n++; } gene[n] = 0; }
}

/* ---------- loading ---------- */

static void load_db(const char *path) {
    size_t n; char *buf = read_all(path, &n);
    DB = (Seq *)malloc(sizeof(Seq) * MAXSEQ);
    char acc[32], gene[32];
    char *p = buf;
    int clen = 0, cap = 0;
    uint8_t *cur = NULL;
    while (p && *p) {
        char *e = strchr(p, '\n');
        if (!e) e = p + strlen(p);
        if (*p == '>') {
            if (clen > 0) {
                DB[NDB].len = clen;
                DB[NDB].s = (uint8_t *)malloc((size_t)clen);
                memcpy(DB[NDB].s, cur, (size_t)clen);
                memcpy(DB[NDB].acc, acc, sizeof(acc));
                memcpy(DB[NDB].gene, gene, sizeof(gene));
                NDB++;
                if (NDB >= MAXSEQ) { fprintf(stderr, "too many seqs\n"); exit(1); }
            }
            char save = *e; *e = 0;
            parse_header(p + 1, acc, gene);
            *e = save;
            clen = 0;
        } else {
            for (char *c = p; c < e; c++) {
                if (*c == '\r' || *c == ' ' || *c == '\t') continue;
                if (clen + 1 > cap) { cap = cap ? cap * 2 : 4096; cur = (uint8_t *)realloc(cur, (size_t)cap); }
                int8_t v = CODE[(unsigned char)*c];
                cur[clen++] = (v < 0) ? (uint8_t)255 : (uint8_t)v;
            }
        }
        p = (*e) ? e + 1 : NULL;
    }
    if (clen > 0) {
        DB[NDB].len = clen;
        DB[NDB].s = (uint8_t *)malloc((size_t)clen);
        memcpy(DB[NDB].s, cur, (size_t)clen);
        memcpy(DB[NDB].acc, acc, sizeof(acc));
        memcpy(DB[NDB].gene, gene, sizeof(gene));
        NDB++;
    }
    free(cur); free(buf);
    fprintf(stderr, "[db] %d sequences\n", NDB);
}

static void load_query(const char *path) {
    size_t n; char *buf = read_all(path, &n);
    char *p = buf; int seen = 0;
    while (p && *p) {
        char *e = strchr(p, '\n');
        if (!e) e = p + strlen(p);
        if (*p == '>') {
            if (!seen) { char save = *e; *e = 0; parse_header(p + 1, QACC, QGENE); *e = save; seen = 1; }
        } else if (seen) {
            for (char *c = p; c < e; c++) {
                if (*c == '\r' || *c == ' ' || *c == '\t') continue;
                int8_t v = CODE[(unsigned char)*c];
                Q[QLEN++] = (v < 0) ? (uint8_t)255 : (uint8_t)v;
            }
        }
        p = (*e) ? e + 1 : NULL;
    }
    free(buf);
    fprintf(stderr, "[query] %s %s len=%d\n", QACC, QGENE, QLEN);
}

/* ---------- composition gate ---------- */

static inline int gate_ok(const uint8_t *w, int L) {
    if (!GATE) return 1;
    int cnt[20]; memset(cnt, 0, sizeof(cnt));
    int distinct = 0;
    for (int i = 0; i < L; i++) {
        int c = w[i];
        if (c >= 20) return 0;
        if (cnt[c]++ == 0) distinct++;
    }
    int mx = 0;
    for (int c = 0; c < 20; c++) if (cnt[c] > mx) mx = cnt[c];
    if (distinct * 100 < 60 * L) return 0;   /* >= 60% distinct residues */
    if (mx * 100 > 40 * L) return 0;         /* no residue > 40% of motif */
    return 1;
}

/* ---------- piece tables ---------- */

typedef struct {
    int      m;
    int      nb;
    int32_t *head;
    int32_t *nxt;
    uint64_t *ent;
    int      n;
} Piece;

static int pow20(int m) { int v = 1; for (int i = 0; i < m; i++) v *= 20; return v; }

typedef struct {
    int    L, d, parts;
    int    sizes[8], offs[8];
    Piece  P[8];
    long   nvalid;          /* number of indexed db windows */
    int    counts[4];       /* counts[mismatch] */
} LTab;

static LTab *TABS[16];
static int   NTAB = 0;

static LTab *build_tab(int L, int d) {
    LTab *t = (LTab *)calloc(1, sizeof(LTab));
    t->L = L; t->d = d; t->parts = d + 1;
    int base = L / t->parts, rem = L % t->parts, off = 0;
    for (int i = 0; i < t->parts; i++) { t->sizes[i] = base + (i < rem ? 1 : 0); t->offs[i] = off; off += t->sizes[i]; }
    /* count valid windows */
    long valid = 0;
    for (int s = 0; s < NDB; s++) {
        int len = DB[s].len;
        if (len < L) continue;
        uint8_t *sq = DB[s].s;
        for (int pos = 0; pos + L <= len; pos++) {
            int ok = 1;
            for (int k = 0; k < L; k++) if (sq[pos + k] >= 20) { ok = 0; break; }
            if (!ok) continue;
            if (GATE) {
                int cnt[20]; memset(cnt, 0, sizeof(cnt));
                int distinct = 0;
                for (int k = 0; k < L; k++) if (cnt[sq[pos + k]]++ == 0) distinct++;
                int mx = 0;
                for (int c = 0; c < 20; c++) if (cnt[c] > mx) mx = cnt[c];
                if (distinct * 100 < 60 * L) continue;
                if (mx * 100 > 40 * L) continue;
            }
            valid++;
        }
    }
    t->nvalid = valid;
    fprintf(stderr, "[tab L=%d d=%d] indexed windows = %ld\n", L, d, valid);
    for (int q = 0; q < t->parts; q++) {
        t->P[q].m = t->sizes[q];
        t->P[q].nb = pow20(t->sizes[q]);
        t->P[q].head = (int32_t *)malloc(sizeof(int32_t) * (size_t)t->P[q].nb);
        for (int i = 0; i < t->P[q].nb; i++) t->P[q].head[i] = -1;
        t->P[q].nxt = (int32_t *)malloc(sizeof(int32_t) * (size_t)(valid + 1));
        t->P[q].ent = (uint64_t *)malloc(sizeof(uint64_t) * (size_t)(valid + 1));
        t->P[q].n = 0;
    }
    for (int s = 0; s < NDB; s++) {
        int len = DB[s].len;
        if (len < L) continue;
        uint8_t *sq = DB[s].s;
        for (int pos = 0; pos + L <= len; pos++) {
            int ok = 1;
            for (int k = 0; k < L; k++) if (sq[pos + k] >= 20) { ok = 0; break; }
            if (!ok) continue;
            if (GATE) {
                int cnt[20]; memset(cnt, 0, sizeof(cnt));
                int distinct = 0;
                for (int k = 0; k < L; k++) if (cnt[sq[pos + k]]++ == 0) distinct++;
                int mx = 0;
                for (int c = 0; c < 20; c++) if (cnt[c] > mx) mx = cnt[c];
                if (distinct * 100 < 60 * L) continue;
                if (mx * 100 > 40 * L) continue;
            }
            uint64_t packed = ((uint64_t)(uint32_t)s << 20) | (uint64_t)(uint32_t)pos;
            for (int q = 0; q < t->parts; q++) {
                int m = t->sizes[q], o = t->offs[q];
                int key = 0;
                for (int k = 0; k < m; k++) key = key * 20 + sq[pos + o + k];
                int idx = t->P[q].n++;
                t->P[q].ent[idx] = packed;
                t->P[q].nxt[idx] = t->P[q].head[key];
                t->P[q].head[key] = idx;
            }
        }
    }
    return t;
}

/* ---------- scan ---------- */

typedef struct {
    int   run;
    int   tabidx;
    const uint8_t *q;
    int   qlen;
    FILE *detail;       /* non-NULL only for run 0 */
    long  nrow;
    long  count[4];
    int   best_nid;
} Task;

static pthread_mutex_t MTX = PTHREAD_MUTEX_INITIALIZER;
static LTab *GLOB_T[16];
static int   NGLOB;

/* per-thread stamp buffers */
static uint32_t *stamp_buf[MAXTHR];
static long      stamp_n[MAXTHR];
static uint32_t  stamp_epoch[MAXTHR];

static void *worker(void *arg) {
    Task *tk = (Task *)arg;
    int tid = tk->tabidx;           /* one slot per table; tabs run concurrently within a run */
    LTab *t = GLOB_T[tk->tabidx];
    int L = t->L, d = t->d;
    const uint8_t *q = tk->q;
    int qlen = tk->qlen;
    long nv = t->nvalid + 1;
    if (stamp_n[tid] < nv) {
        free(stamp_buf[tid]);
        stamp_buf[tid] = (uint32_t *)calloc((size_t)nv, sizeof(uint32_t));
        stamp_n[tid] = nv;
        stamp_epoch[tid] = 0;
    }
    uint32_t *stamp = stamp_buf[tid];
    uint64_t epoch = ++stamp_epoch[tid];
    if (epoch > 4000000000ULL) { memset(stamp, 0, sizeof(uint32_t) * (size_t)stamp_n[tid]); epoch = stamp_epoch[tid] = 1; }

    long cnt[4] = {0, 0, 0, 0};
    long rows = 0;

    for (int i = 0; i + L <= qlen; i++) {
        const uint8_t *w = q + i;
        if (!gate_ok(w, L)) continue;
        for (int p = 0; p < t->parts; p++) {
            int m = t->sizes[p], o = t->offs[p];
            int key = 0, bad = 0;
            for (int k = 0; k < m; k++) { int c = w[o + k]; if (c >= 20) { bad = 1; break; } key = key * 20 + c; }
            if (bad) continue;
            for (int32_t e = t->P[p].head[key]; e >= 0; e = t->P[p].nxt[e]) {
                if (stamp[e] == epoch) continue;
                uint64_t packed = t->P[p].ent[e];
                int s = (int)(packed >> 20);
                int pos = (int)(packed & 0xFFFFF);
                uint8_t *ts = DB[s].s + pos;
                int nid = 0, mis = 0;
                for (int k = 0; k < L; k++) {
                    if (w[k] == ts[k]) nid++;
                    else if (++mis > d) break;
                }
                if (mis > d) continue;
                stamp[e] = epoch;
                cnt[L - nid] ++;
                rows++;
                if (tk->detail && rows <= ROWCAP) {
                    char bq[16], bt[16];
                    for (int k = 0; k < L; k++) { bq[k] = AALET[w[k]]; bt[k] = AALET[ts[k]]; }
                    bq[L] = 0; bt[L] = 0;
                    fprintf(tk->detail, "%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n",
                            L, i + 1, DB[s].acc, DB[s].gene, pos + 1, nid, bq, bt);
                }
            }
        }
    }
    tk->nrow = rows;
    for (int k = 0; k <= d; k++) tk->count[k] = cnt[k];
    return NULL;
}

/* ---------- shuffles ---------- */

static void shuffle_windows(uint8_t *dst, const uint8_t *src, int len, int w, uint64_t seed) {
    memcpy(dst, src, (size_t)len);
    uint64_t st = seed ? seed : 88172645463325252ULL;
    for (int a = 0; a < len; a += w) {
        int b = a + w; if (b > len) b = len;
        for (int i = b - 1; i > a; i--) {
            st ^= st << 13; st ^= st >> 7; st ^= st << 17;
            int j = a + (int)(st % (uint64_t)(i - a + 1));
            uint8_t tmp = dst[i]; dst[i] = dst[j]; dst[j] = tmp;
        }
    }
}

/* ---------- main ---------- */

int main(int argc, char **argv) {
    const char *qpath = NULL, *dpath = NULL, *prefix = "motif";
    int null_b = 50, seed = 42, win = 20;
    int Ls[16], nL = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--query") && i + 1 < argc) qpath = argv[++i];
        else if (!strcmp(argv[i], "--db") && i + 1 < argc) dpath = argv[++i];
        else if (!strcmp(argv[i], "--out-prefix") && i + 1 < argc) prefix = argv[++i];
        else if (!strcmp(argv[i], "--null-b") && i + 1 < argc) null_b = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--gate") && i + 1 < argc) GATE = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--window") && i + 1 < argc) win = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--threads") && i + 1 < argc) THREADS = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--lens") && i + 1 < argc) {
            char *s = argv[++i];
            char *tok = strtok(s, ",");
            while (tok) { Ls[nL++] = atoi(tok); tok = strtok(NULL, ","); }
        } else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 2; }
    }
    if (!qpath || !dpath) { fprintf(stderr, "need --query and --db\n"); return 2; }
    if (nL == 0) { Ls[nL++] = 6; Ls[nL++] = 7; Ls[nL++] = 8; Ls[nL++] = 9; Ls[nL++] = 10; }
    if (THREADS > MAXTHR) THREADS = MAXTHR;

    init_code();
    load_db(dpath);
    load_query(qpath);

    fprintf(stderr, "[cfg] gate=%d null_b=%d seed=%d window=%d threads=%d\n", GATE, null_b, seed, win, THREADS);

    /* build tables */
    NTAB = 0;
    for (int i = 0; i < nL; i++) {
        int L = Ls[i];
        int d = (L <= 6) ? 1 : ((L <= 9) ? 2 : 3);
        GLOB_T[NTAB] = build_tab(L, d);
        NTAB++;
    }

    /* prepare queries: run 0 = real, runs 1..null_b = window shuffle */
    int NRUN = null_b + 1;
    uint8_t **QRY = (uint8_t **)malloc(sizeof(uint8_t *) * (size_t)NRUN);
    QRY[0] = (uint8_t *)malloc((size_t)QLEN);
    memcpy(QRY[0], Q, (size_t)QLEN);
    for (int r = 1; r < NRUN; r++) {
        QRY[r] = (uint8_t *)malloc((size_t)QLEN);
        shuffle_windows(QRY[r], Q, QLEN, win, (uint64_t)seed * 1000003ULL + (uint64_t)r);
    }

    /* counts[r][tab][tier] */
    long ***counts = (long ***)malloc(sizeof(long **) * (size_t)NRUN);
    for (int r = 0; r < NRUN; r++) {
        counts[r] = (long **)malloc(sizeof(long *) * (size_t)NTAB);
        for (int t = 0; t < NTAB; t++) {
            counts[r][t] = (long *)calloc(4, sizeof(long));
        }
    }

    /* run 0 detail files per table */
    FILE **det = (FILE **)calloc((size_t)NTAB, sizeof(FILE *));
    char fn[512];

    /* execute run 0 first (parallel over tables) */
    {
        pthread_t th[16]; Task *tk = (Task *)calloc((size_t)NTAB, sizeof(Task));
        for (int t = 0; t < NTAB; t++) {
            snprintf(fn, sizeof(fn), "%s_real_L%d.tsv", prefix, GLOB_T[t]->L);
            det[t] = fopen(fn, "wb");
            fprintf(det[t], "L\tqstart\tacc\tgene\ttstart\tn_ident\tmotif_q\tmotif_t\n");
            tk[t].run = 0; tk[t].tabidx = t; tk[t].q = QRY[0]; tk[t].qlen = QLEN;
            tk[t].detail = det[t];
        }
        for (int t = 0; t < NTAB; t++) pthread_create(&th[t], NULL, worker, &tk[t]);
        for (int t = 0; t < NTAB; t++) pthread_join(th[t], NULL);
        for (int t = 0; t < NTAB; t++) {
            fclose(det[t]);
            for (int k = 0; k < 4; k++) counts[0][t][k] = tk[t].count[k];
            fprintf(stderr, "[real L=%d d=%d] rows=%ld  mm0=%ld mm1=%ld mm2=%ld mm3=%ld\n", GLOB_T[t]->L, GLOB_T[t]->d,
                    tk[t].nrow, tk[t].count[0], tk[t].count[1], tk[t].count[2], tk[t].count[3]);
        }
        free(tk);
    }

    /* null runs */
    for (int r = 1; r < NRUN; r++) {
        pthread_t th[16]; Task *tk = (Task *)calloc((size_t)NTAB, sizeof(Task));
        for (int t = 0; t < NTAB; t++) {
            tk[t].run = r; tk[t].tabidx = t; tk[t].q = QRY[r]; tk[t].qlen = QLEN; tk[t].detail = NULL;
        }
        for (int t = 0; t < NTAB; t++) pthread_create(&th[t], NULL, worker, &tk[t]);
        for (int t = 0; t < NTAB; t++) pthread_join(th[t], NULL);
        for (int t = 0; t < NTAB; t++) for (int k = 0; k < 4; k++) counts[r][t][k] = tk[t].count[k];
        free(tk);
        if ((r % 10) == 0) fprintf(stderr, "[null] %d/%d done\n", r, null_b);
    }

    /* dump counts */
    snprintf(fn, sizeof(fn), "%s_counts.tsv", prefix);
    FILE *cf = fopen(fn, "wb");
    fprintf(cf, "run\tL\tmismatch\tcount\n");
    for (int r = 0; r < NRUN; r++)
        for (int t = 0; t < NTAB; t++)
            for (int k = 0; k < 4; k++)
                if (counts[r][t][k] > 0)
                    fprintf(cf, "%d\t%d\t%d\t%ld\n", r, GLOB_T[t]->L, k, counts[r][t][k]);
    fclose(cf);
    fprintf(stderr, "[done] wrote %s\n", fn);
    return 0;
}
