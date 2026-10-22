/* ===========================================================================
 * mimic_scan.c — faithful C port of the Biomni ebna1-mimicscan v1.3 core.
 *
 * Ported 1:1 from:
 *   R/seeds.R        build_kmer_index() / kmer_selfscore() / gar_fraction()
 *                    / select_seeds()
 *   R/scan.R         scan_proteome()  (dedup: by acc,qstart,qend then
 *                    acc,tstart,tend, after ordering by -gscore)
 *                    scan_proteome_exhaustive()
 *   R/scan_core.cpp  scan_seeds() / sw_one_masked() / extend_to_min_len()
 *
 * Deliberately preserved (so the reproduction can be compared number-for-number
 * with the R package):
 *   - BLOSUM62 24x24 in aaidx order A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,B,Z,X,*
 *   - k = 4, seed self-score >= 14, GAR fraction < 0.75
 *   - ungapped X-drop extension (X = 6) around every seed occurrence
 *   - gapped Smith-Waterman inside [seed +/- slack] window (slack = 20),
 *     LINEAR gap cost = gap_open per residue (gap_ext is never used by the
 *     original code — this is a known dead parameter, reproduced on purpose)
 *   - strict SW (H never < 0), tie order diagonal > up > left
 *   - min_len diagonal "rescue" that extends a too-short optimum, keeping the
 *     pre-extension core coordinates as separate columns
 *   - exhaustive seedless mode: full SW per protein, query-footprint masking,
 *     up to n_iter distinct alignments per protein
 *
 * Differences from the original (documented, not silent):
 *   - the k-mer index is a counting-sort bucket table instead of a keyed
 *     data.table; selection semantics are identical
 *   - sort/dedup tie-breaking is deterministic but not byte-identical to R's
 *     data.table ordering
 *
 * Build:  gcc -O3 -march=native -pthread -o mimic_scan mimic_scan.c
 * Usage:  mimic_scan --query Q.fasta --db DB.fasta --mode seeded|exhaustive
 *                    --out OUT.tsv [params...]
 * =========================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <sys/types.h>
#include <pthread.h>

#define NA_INT INT_MIN

/* --------------------------------------------------------------------------
 * BLOSUM62 (values transcribed from R/blosum62.R; row/col 23 = '*' filled
 * the same way the R code fills it)
 * -------------------------------------------------------------------------- */
static int S[24][24];

static const char AA24[] = "ARNDCQEGHILKMFPSTWYVBZX*";

static const int BL62[23][23] = {
/*          A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X */
/*A*/ {     4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0},
/*R*/ {    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1},
/*N*/ {    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1},
/*D*/ {    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1},
/*C*/ {     0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2},
/*Q*/ {    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1},
/*E*/ {    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1},
/*G*/ {     0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1},
/*H*/ {    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1},
/*I*/ {    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1},
/*L*/ {    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1},
/*K*/ {    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1},
/*M*/ {    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1},
/*F*/ {    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1},
/*P*/ {    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2},
/*S*/ {     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0},
/*T*/ {     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0},
/*W*/ {    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2},
/*Y*/ {    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1},
/*V*/ {     0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1},
/*B*/ {    -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1},
/*Z*/ {    -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1},
/*X*/ {     0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1}
};

static void init_matrix(void) {
  for (int i = 0; i < 23; i++)
    for (int j = 0; j < 23; j++)
      S[i][j] = BL62[i][j];
  for (int k = 0; k < 24; k++) { S[23][k] = -4; S[k][23] = -4; }
  S[23][23] = 1;
}

static inline int aaidx(char c) {
  switch (c) {
    case 'A': return 0;  case 'R': return 1;  case 'N': return 2;  case 'D': return 3;
    case 'C': return 4;  case 'Q': return 5;  case 'E': return 6;  case 'G': return 7;
    case 'H': return 8;  case 'I': return 9;  case 'L': return 10; case 'K': return 11;
    case 'M': return 12; case 'F': return 13; case 'P': return 14; case 'S': return 15;
    case 'T': return 16; case 'W': return 17; case 'Y': return 18; case 'V': return 19;
    case 'B': return 20; case 'Z': return 21; case 'X': return 22;
    case 'U': return 22; case 'O': return 22; default: return 22;
  }
}

/* --------------------------------------------------------------------------
 * FASTA reading
 * -------------------------------------------------------------------------- */
typedef struct {
  char *acc;      /* UniProt accession (field 2 of sp|ACC|NAME)           */
  char *gene;     /* GN=... if present, else entry name                    */
  char *seq;      /* raw sequence (uppercase, no whitespace)              */
  unsigned char *enc;  /* aaidx codes                                      */
  int len;
} Rec;

typedef struct { Rec *v; int n, cap; } Recs;

static void recs_push(Recs *r, Rec rec) {
  if (r->n == r->cap) {
    r->cap = r->cap ? r->cap * 2 : 1024;
    r->v = (Rec *)realloc(r->v, (size_t)r->cap * sizeof(Rec));
    if (!r->v) { fprintf(stderr, "OOM\n"); exit(1); }
  }
  r->v[r->n++] = rec;
}

static char *xstrdup(const char *s) {
  size_t n = strlen(s) + 1;
  char *p = (char *)malloc(n);
  if (!p) { fprintf(stderr, "OOM\n"); exit(1); }
  memcpy(p, s, n);
  return p;
}

/* header: >sp|P02452|CO1A1_HUMAN Description ... GN=COL1A1 PE=1 SV=5 */
static void parse_header(const char *hdr, char **acc, char **gene) {
  char accbuf[128] = {0};
  char genebuf[128] = {0};
  /* accession: text between first and second '|', else first token */
  const char *p1 = strchr(hdr, '|');
  if (p1) {
    const char *p2 = strchr(p1 + 1, '|');
    if (p2 && (size_t)(p2 - p1 - 1) < sizeof(accbuf)) {
      memcpy(accbuf, p1 + 1, (size_t)(p2 - p1 - 1));
    }
  }
  if (!accbuf[0]) {
    size_t k = 0;
    while (hdr[k] && hdr[k] != ' ' && k < sizeof(accbuf) - 1) { accbuf[k] = hdr[k]; k++; }
    accbuf[k] = 0;
  }
  /* gene: GN= token */
  const char *g = strstr(hdr, "GN=");
  if (g) {
    g += 3;
    size_t k = 0;
    while (g[k] && g[k] != ' ' && g[k] != ';' && k < sizeof(genebuf) - 1) { genebuf[k] = g[k]; k++; }
    genebuf[k] = 0;
  }
  if (!genebuf[0]) {
    /* entry name = field 3 */
    if (p1) {
      const char *p2 = strchr(p1 + 1, '|');
      if (p2) {
        const char *p3 = p2 + 1;
        size_t k = 0;
        while (p3[k] && p3[k] != ' ' && k < sizeof(genebuf) - 1) { genebuf[k] = p3[k]; k++; }
        genebuf[k] = 0;
      }
    }
  }
  *acc = xstrdup(accbuf);
  *gene = xstrdup(genebuf[0] ? genebuf : accbuf);
}

static void read_fasta(const char *path, Recs *out) {
  FILE *f = fopen(path, "r");
  if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
  char *line = NULL; size_t cap = 0; ssize_t n;
  char *acc = NULL, *gene = NULL, *seq = NULL;
  size_t slen = 0, scap = 0;
  while ((n = getline(&line, &cap, f)) != -1) {
    while (n > 0 && (line[n-1] == '\n' || line[n-1] == '\r')) line[--n] = 0;
    if (n == 0) continue;
    if (line[0] == '>') {
      if (seq) {
        Rec r; r.acc = acc; r.gene = gene; r.seq = seq; r.len = (int)slen;
        r.enc = (unsigned char *)malloc((size_t)slen + 1);
        for (size_t i = 0; i < slen; i++) r.enc[i] = (unsigned char)aaidx(seq[i]);
        recs_push(out, r);
      }
      parse_header(line + 1, &acc, &gene);
      seq = NULL; slen = 0; scap = 0;
    } else {
      if (slen + (size_t)n + 1 > scap) {
        scap = (slen + (size_t)n + 1) * 2;
        seq = (char *)realloc(seq, scap);
        if (!seq) { fprintf(stderr, "OOM\n"); exit(1); }
      }
      for (ssize_t i = 0; i < n; i++) {
        char c = line[i];
        if (c >= 'a' && c <= 'z') c = (char)(c - 32);
        seq[slen++] = c;
      }
      seq[slen] = 0;
    }
  }
  if (seq) {
    Rec r; r.acc = acc; r.gene = gene; r.seq = seq; r.len = (int)slen;
    r.enc = (unsigned char *)malloc((size_t)slen + 1);
    for (size_t i = 0; i < slen; i++) r.enc[i] = (unsigned char)aaidx(seq[i]);
    recs_push(out, r);
  }
  free(line);
  fclose(f);
}

/* --------------------------------------------------------------------------
 * Result rows
 * -------------------------------------------------------------------------- */
typedef struct {
  int mode;           /* 0 seeded, 1 exhaustive */
  int pi;             /* 1-based protein index */
  int qs, qe, ts, te;
  int alen, nid;
  int uscore, gscore;
  int cqs, cqe, cts, cte, cscore, cnid, clen;
  char ext[32];
  char *gq, *gt;      /* gapped alignment strings (may be empty) */
} Row;

typedef struct { Row *v; int n, cap; } Rows;

static void rows_push(Rows *r, const Row *row) {
  if (r->n == r->cap) {
    r->cap = r->cap ? r->cap * 2 : 4096;
    r->v = (Row *)realloc(r->v, (size_t)r->cap * sizeof(Row));
    if (!r->v) { fprintf(stderr, "OOM\n"); exit(1); }
  }
  r->v[r->n++] = *row;
}

/* --------------------------------------------------------------------------
 * extend_to_min_len — verbatim port of scan_core.cpp
 * -------------------------------------------------------------------------- */
static int extend_to_min_len(const char *query, const char *subj,
                             int qlo, int slo2, int nq, int ns,
                             int *i, int *j, int *gqi, int *gsj,
                             char *gqa, char *gsa, int *gqa_len, int *gsa_len,
                             int *gnid, int *galen, int *curscore,
                             int min_score, int min_len) {
  while (*galen < min_len) {
    int addL = -1000000, addR = -1000000;
    if (*i > 0 && *j > 0)
      addL = S[aaidx(query[qlo + *i - 1])][aaidx(subj[slo2 + *j - 1])];
    if (*gqi < nq && *gsj < ns)
      addR = S[aaidx(query[qlo + *gqi])][aaidx(subj[slo2 + *gsj])];
    int useL = (addL >= addR);
    int add = useL ? addL : addR;
    if (*curscore + add < min_score) return 0;
    if (useL) {
      char a = query[qlo + *i - 1], b = subj[slo2 + *j - 1];
      memmove(gqa + 1, gqa, (size_t)*gqa_len + 1); gqa[0] = a; (*gqa_len)++;
      memmove(gsa + 1, gsa, (size_t)*gsa_len + 1); gsa[0] = b; (*gsa_len)++;
      if (a == b) (*gnid)++;
      (*i)--; (*j)--;
    } else {
      char a = query[qlo + *gqi], b = subj[slo2 + *gsj];
      gqa[(*gqa_len)++] = a; gqa[*gqa_len] = 0;
      gsa[(*gsa_len)++] = b; gsa[*gsa_len] = 0;
      if (a == b) (*gnid)++;
      (*gqi)++; (*gsj)++;
    }
    (*galen)++; *curscore += add;
  }
  return 1;
}

/* --------------------------------------------------------------------------
 * Parameters
 * -------------------------------------------------------------------------- */
typedef struct {
  int kmer, seed_score_min;
  double seed_gar_max;
  int xdrop, slack, gap_open, gap_ext, min_len;
  int min_score;
  int n_iter, threads;
} Params;

/* --------------------------------------------------------------------------
 * k-mer index (counting sort over 24^k buckets, k = 4)
 * -------------------------------------------------------------------------- */
typedef struct {
  int k;
  uint32_t nbuckets;
  uint32_t *start;   /* nbuckets + 1 */
  uint64_t *items;   /* (pi << 32) | pos   (both 0-based) */
  uint32_t nitems;
} KIdx;

static inline uint32_t kcode(const unsigned char *e, int i, int k) {
  uint32_t c = 0;
  for (int t = 0; t < k; t++) c = c * 24u + (uint32_t)e[i + t];
  return c;
}

static void build_index(KIdx *idx, const Recs *prot, int k) {
  uint32_t nb = 1;
  for (int t = 0; t < k; t++) nb *= 24u;
  idx->k = k; idx->nbuckets = nb;
  uint32_t *cnt = (uint32_t *)calloc((size_t)nb + 1, sizeof(uint32_t));
  if (!cnt) { fprintf(stderr, "OOM\n"); exit(1); }
  uint64_t total = 0;
  for (int p = 0; p < prot->n; p++) {
    const Rec *r = &prot->v[p];
    if (r->len < k) continue;
    for (int i = 0; i + k <= r->len; i++) cnt[kcode(r->enc, i, k)]++;
    total += (uint64_t)(r->len - k + 1);
  }
  uint32_t *start = (uint32_t *)malloc(((size_t)nb + 1) * sizeof(uint32_t));
  uint64_t acc = 0;
  for (uint32_t b = 0; b < nb; b++) { start[b] = (uint32_t)acc; acc += cnt[b]; }
  start[nb] = (uint32_t)acc;
  memcpy(cnt, start, (size_t)nb * sizeof(uint32_t));   /* cursor */
  uint64_t *items = (uint64_t *)malloc((size_t)total * sizeof(uint64_t));
  if (!items) { fprintf(stderr, "OOM\n"); exit(1); }
  for (int p = 0; p < prot->n; p++) {
    const Rec *r = &prot->v[p];
    if (r->len < k) continue;
    for (int i = 0; i + k <= r->len; i++) {
      uint32_t c = kcode(r->enc, i, k);
      items[cnt[c]++] = ((uint64_t)(uint32_t)p << 32) | (uint32_t)i;
    }
  }
  free(cnt);
  idx->start = start;
  idx->items = items;
  idx->nitems = (uint32_t)total;
}

/* --------------------------------------------------------------------------
 * Seeded scan (port of scan_seeds)
 * -------------------------------------------------------------------------- */
static void scan_seeded(const Recs *prot, const Rec *query, const KIdx *idx,
                        const Params *P, Rows *out) {
  const int k = P->kmer;
  const unsigned char *qenc = query->enc;
  const int qlen = query->len;
  if (qlen < k) return;

  int nseeds = 0;
  for (int e = 0; e + k <= qlen; e++) {
    /* kmer_selfscore: sum of diagonal */
    int ss = 0;
    for (int t = 0; t < k; t++) { int a = qenc[e + t]; ss += S[a][a]; }
    if (ss < P->seed_score_min) continue;
    /* gar_fraction */
    int gar = 0;
    for (int t = 0; t < k; t++) { int a = qenc[e + t]; if (a == 7 || a == 0 || a == 1) gar++; }
    double gf = (double)gar / (double)k;
    if (!(gf < P->seed_gar_max)) continue;
    nseeds++;
  }
  fprintf(stderr, "[seed] %d of %d query k-mers pass self-score>=%d & GAR<%.2f\n",
          nseeds, qlen - k + 1, P->seed_score_min, P->seed_gar_max);

  /* working buffers for one window alignment */
  enum { WBUF = 2048 };
  int Hbuf[64][64];
  char Pbuf[64][64];
  char gqa[WBUF], gsa[WBUF];

  long nocc = 0;
  for (int e = 0; e + k <= qlen; e++) {
    int ss = 0;
    for (int t = 0; t < k; t++) { int a = qenc[e + t]; ss += S[a][a]; }
    if (ss < P->seed_score_min) continue;
    int gar = 0;
    for (int t = 0; t < k; t++) { int a = qenc[e + t]; if (a == 7 || a == 0 || a == 1) gar++; }
    if (!((double)gar / (double)k < P->seed_gar_max)) continue;

    uint32_t c = kcode(qenc, e, k);
    for (uint32_t s = idx->start[c]; s < idx->start[c + 1]; s++) {
      int pi = (int)(idx->items[s] >> 32);
      int spos = (int)(idx->items[s] & 0xffffffffu);
      const Rec *subj = &prot->v[pi];
      const int slen = subj->len;
      if (spos + k > slen) continue;
      nocc++;

      const unsigned char *senc = subj->enc;

      /* ---------- ungapped X-drop extension ---------- */
      int lo = e, hi = e + k - 1;
      int slo = spos, shi = spos + k - 1;
      int seed_score = 0;
      for (int t = 0; t < k; t++) seed_score += S[qenc[e + t]][senc[spos + t]];
      int best = seed_score, blo = lo, bhi = hi, bslo = slo, bshi = shi;
      int cur = best;
      for (int i = lo - 1, j = slo - 1; i >= 0 && j >= 0; i--, j--) {
        cur += S[qenc[i]][senc[j]];
        if (cur > best) { best = cur; blo = i; bslo = j; }
        if (cur < best - P->xdrop) break;
      }
      cur = best;
      for (int i = bhi + 1, j = bshi + 1; i < qlen && j < slen; i++, j++) {
        cur += S[qenc[i]][senc[j]];
        if (cur > best) { best = cur; bhi = i; bshi = j; }
        if (cur < best - P->xdrop) break;
      }
      int ulen = bhi - blo + 1;
      if (ulen >= P->min_len) {
        int nid = 0;
        for (int i = blo; i <= bhi; i++)
          if (query->seq[i] == subj->seq[bslo + (i - blo)]) nid++;
        Row r; memset(&r, 0, sizeof r);
        r.mode = 0; r.pi = pi + 1;
        r.qs = blo + 1; r.qe = bhi + 1; r.ts = bslo + 1; r.te = bshi + 1;
        r.alen = ulen; r.nid = nid; r.uscore = best; r.gscore = NA_INT;
        r.cqs = r.qs; r.cqe = r.qe; r.cts = r.ts; r.cte = r.te;
        r.cscore = best; r.cnid = nid; r.clen = ulen;
        strcpy(r.ext, "none");
        r.gq = xstrdup(""); r.gt = xstrdup("");
        rows_push(out, &r);
      }

      /* ---------- gapped SW in [seed +/- slack] window ---------- */
      int qlo = e - P->slack; if (qlo < 0) qlo = 0;
      int qhi = e + k - 1 + P->slack; if (qhi > qlen - 1) qhi = qlen - 1;
      int slo2 = spos - P->slack; if (slo2 < 0) slo2 = 0;
      int shi2 = spos + k - 1 + P->slack; if (shi2 > slen - 1) shi2 = slen - 1;
      int nq = qhi - qlo + 1, ns = shi2 - slo2 + 1;
      if (nq < 3 || ns < 3) continue;
      if (nq >= 64 || ns >= 64) { fprintf(stderr, "window too large (%d,%d)\n", nq, ns); continue; }

      int gmax = 0, gqi = 0, gsj = 0;
      for (int i = 1; i <= nq; i++) {
        Hbuf[i][0] = 0; Pbuf[i][0] = 0;
      }
      for (int j = 0; j <= ns; j++) { Hbuf[0][j] = 0; Pbuf[0][j] = 0; }
      for (int i = 1; i <= nq; i++) {
        int ai = qenc[qlo + i - 1];
        for (int j = 1; j <= ns; j++) {
          int sc = S[ai][senc[slo2 + j - 1]];
          int d = Hbuf[i-1][j-1] + sc;
          int u = Hbuf[i-1][j] - P->gap_open;
          int l = Hbuf[i][j-1] - P->gap_open;
          int v = 0; char p = 0;
          if (d >= u && d >= l && d >= v) { v = d; p = 1; }
          else if (u >= l && u >= v)      { v = u; p = 2; }
          else if (l >= v)                { v = l; p = 3; }
          Hbuf[i][j] = v; Pbuf[i][j] = p;
          if (v > gmax) { gmax = v; gqi = i; gsj = j; }
        }
      }
      if (gmax < P->min_score) continue;

      /* traceback */
      int gqa_len = 0, gsa_len = 0, gnid = 0, galen = 0;
      gqa[0] = 0; gsa[0] = 0;
      int i = gqi, j = gsj;
      /* build reversed then flip */
      char rq[WBUF], rs[WBUF];
      int m = 0;
      while (i > 0 && j > 0 && Pbuf[i][j] != 0) {
        if (Pbuf[i][j] == 1) {
          char a = query->seq[qlo + i - 1], b = subj->seq[slo2 + j - 1];
          rq[m] = a; rs[m] = b; m++;
          if (a == b) gnid++;
          galen++; i--; j--;
        } else if (Pbuf[i][j] == 2) {
          rq[m] = query->seq[qlo + i - 1]; rs[m] = '-'; m++;
          galen++; i--;
        } else {
          rq[m] = '-'; rs[m] = subj->seq[slo2 + j - 1]; m++;
          galen++; j--;
        }
        if (m >= WBUF - 2) break;
      }
      for (int t = 0; t < m; t++) { gqa[t] = rq[m - 1 - t]; gsa[t] = rs[m - 1 - t]; }
      gqa[m] = 0; gsa[m] = 0;
      gqa_len = m; gsa_len = m;

      int c_qs, c_qe, c_ss, c_se, c_score, c_nid, c_len;
      char extm[32] = "none";
      c_qs = qlo + i + 1; c_qe = qlo + gqi; c_ss = slo2 + j + 1; c_se = slo2 + gsj;
      c_score = gmax; c_nid = gnid; c_len = galen;
      if (galen < P->min_len) {
        int ok = extend_to_min_len(query->seq, subj->seq, qlo, slo2, nq, ns,
                                   &i, &j, &gqi, &gsj, gqa, gsa, &gqa_len, &gsa_len,
                                   &gnid, &galen, &gmax, P->min_score, P->min_len);
        if (!ok) continue;
        strcpy(extm, "diagonal_min_len_rescue");
        gqa_len = (int)strlen(gqa); gsa_len = (int)strlen(gsa);
      }

      Row r; memset(&r, 0, sizeof r);
      r.mode = 0; r.pi = pi + 1;
      r.qs = qlo + i + 1; r.qe = qlo + gqi;
      r.ts = slo2 + j + 1; r.te = slo2 + gsj;
      r.alen = galen; r.nid = gnid;
      r.uscore = NA_INT; r.gscore = gmax;
      r.cqs = c_qs; r.cqe = c_qe; r.cts = c_ss; r.cte = c_se;
      r.cscore = c_score; r.cnid = c_nid; r.clen = c_len;
      strcpy(r.ext, extm);
      r.gq = xstrdup(gqa); r.gt = xstrdup(gsa);
      rows_push(out, &r);
    }
  }
  fprintf(stderr, "[seed] %ld seed occurrences extended\n", nocc);
}

/* --------------------------------------------------------------------------
 * Exhaustive scan (port of sw_one_masked / scan_exhaustive)
 * -------------------------------------------------------------------------- */
typedef struct {
  const Recs *prot;
  const Rec *query;
  const Params *P;
  int tid, nthreads;
  Rows rows;
} ExhJob;

static void sw_one_masked(const Rec *query, const Rec *subj, int pi1,
                          const Params *P, const unsigned char *qmask,
                          Row *out, int *H, char *Pt) {
  const int qlen = query->len, slen = subj->len;
  out->pi = -1;
  if (slen < 3) return;
  const int nscol = slen + 1;
  memset(H, 0, (size_t)(qlen + 1) * (size_t)nscol * sizeof(int));
  memset(Pt, 0, (size_t)(qlen + 1) * (size_t)nscol);
  int gmax = 0, gqi = 0, gsj = 0;
  for (int i = 1; i <= qlen; i++) {
    if (qmask[i - 1]) continue;
    int ai = query->enc[i - 1];
    int *Hrow = H + (size_t)i * nscol;
    int *Hprev = H + (size_t)(i - 1) * nscol;
    char *Prow = Pt + (size_t)i * nscol;
    for (int j = 1; j <= slen; j++) {
      int sc = S[ai][subj->enc[j - 1]];
      int d = Hprev[j - 1] + sc;
      int u = Hprev[j] - P->gap_open;
      int l = Hrow[j - 1] - P->gap_open;
      int v = 0; char p = 0;
      if (d >= u && d >= l && d >= 0) { v = d; p = 1; }
      else if (u >= l && u >= 0)      { v = u; p = 2; }
      else if (l >= 0)                { v = l; p = 3; }
      Hrow[j] = v; Prow[j] = p;
      if (v > gmax) { gmax = v; gqi = i; gsj = j; }
    }
  }
  if (gmax < P->min_score) return;

  int m = 0, gnid = 0, galen = 0;
  char rq[4096], rs[4096];
  int i = gqi, j = gsj;
  while (i > 0 && j > 0 && Pt[(size_t)i * nscol + j] != 0) {
    char p = Pt[(size_t)i * nscol + j];
    if (p == 1) {
      char a = query->seq[i - 1], b = subj->seq[j - 1];
      rq[m] = a; rs[m] = b; m++;
      if (a == b) gnid++;
      galen++; i--; j--;
    } else if (p == 2) {
      rq[m] = query->seq[i - 1]; rs[m] = '-'; m++;
      galen++; i--;
    } else {
      rq[m] = '-'; rs[m] = subj->seq[j - 1]; m++;
      galen++; j--;
    }
    if (m >= 4000) break;
  }
  char gqa[4096], gsa[4096];
  for (int t = 0; t < m; t++) { gqa[t] = rq[m - 1 - t]; gsa[t] = rs[m - 1 - t]; }
  gqa[m] = 0; gsa[m] = 0;

  int c_qs = i + 1, c_qe = gqi, c_ss = j + 1, c_se = gsj;
  int c_score = gmax, c_nid = gnid, c_len = galen;
  char extm[32] = "none";
  if (galen < P->min_len) {
    int gqa_len = m, gsa_len = m;
    int ok = extend_to_min_len(query->seq, subj->seq, 0, 0, qlen, slen,
                               &i, &j, &gqi, &gsj, gqa, gsa, &gqa_len, &gsa_len,
                               &gnid, &galen, &gmax, P->min_score, P->min_len);
    if (!ok) return;
    strcpy(extm, "diagonal_min_len_rescue");
  }

  memset(out, 0, sizeof *out);
  out->mode = 1; out->pi = pi1;
  out->qs = i + 1; out->qe = gqi; out->ts = j + 1; out->te = gsj;
  out->alen = galen; out->nid = gnid;
  out->uscore = NA_INT; out->gscore = gmax;
  out->cqs = c_qs; out->cqe = c_qe; out->cts = c_ss; out->cte = c_se;
  out->cscore = c_score; out->cnid = c_nid; out->clen = c_len;
  strcpy(out->ext, extm);
  out->gq = xstrdup(gqa); out->gt = xstrdup(gsa);
}

static void *exh_worker(void *arg) {
  ExhJob *J = (ExhJob *)arg;
  const Recs *prot = J->prot;
  const Rec *q = J->query;
  const Params *P = J->P;
  int maxlen = 0;
  for (int p = 0; p < prot->n; p++) if (prot->v[p].len > maxlen) maxlen = prot->v[p].len;
  int nscol = maxlen + 1;
  int *H = (int *)malloc((size_t)(q->len + 1) * (size_t)nscol * sizeof(int));
  char *Pt = (char *)malloc((size_t)(q->len + 1) * (size_t)nscol);
  if (!H || !Pt) { fprintf(stderr, "OOM (exhaustive DP)\n"); exit(1); }
  unsigned char *qmask = (unsigned char *)malloc((size_t)q->len + 1);

  for (int p = J->tid; p < prot->n; p += J->nthreads) {
    memset(qmask, 0, (size_t)q->len);
    for (int it = 0; it < P->n_iter; it++) {
      Row r;
      sw_one_masked(q, &prot->v[p], p + 1, P, qmask, &r, H, Pt);
      if (r.pi <= 0) break;
      rows_push(&J->rows, &r);
      for (int t = r.qs - 1; t <= r.qe - 1; t++) {
        if (t >= 0 && t < q->len) qmask[t] = 1;
      }
    }
  }
  free(H); free(Pt); free(qmask);
  return NULL;
}

/* --------------------------------------------------------------------------
 * Dedup (port of scan_proteome's two unique() passes)
 * -------------------------------------------------------------------------- */
static int cmp_score(const void *a, const void *b) {
  const Row *x = (const Row *)a, *y = (const Row *)b;
  int xg = (x->gscore != NA_INT), yg = (y->gscore != NA_INT);
  if (xg != yg) return yg - xg;               /* gapped rows first */
  if (xg) {
    if (x->gscore != y->gscore) return (y->gscore > x->gscore) ? 1 : -1;
  } else {
    if (x->uscore != y->uscore) return (y->uscore > x->uscore) ? 1 : -1;
  }
  return 0;
}

typedef struct { uint64_t *tab; size_t mask; } HSet;

static void hset_init(HSet *h, size_t n) {
  size_t cap = 1; while (cap < n * 2) cap <<= 1;
  h->tab = (uint64_t *)calloc(cap, sizeof(uint64_t));
  h->mask = cap - 1;
}
static inline int hset_add(HSet *h, uint64_t key) {
  uint64_t k = key ? key : 1;
  size_t i = (size_t)((k * 11400714819323198485ull) >> 40) & h->mask;
  while (h->tab[i]) {
    if (h->tab[i] == k) return 0;
    i = (i + 1) & h->mask;
  }
  h->tab[i] = k;
  return 1;
}

static int dedup(Rows *r) {
  qsort(r->v, (size_t)r->n, sizeof(Row), cmp_score);
  HSet h1, h2;
  hset_init(&h1, (size_t)r->n + 16);
  hset_init(&h2, (size_t)r->n + 16);
  int keep = 0;
  for (int i = 0; i < r->n; i++) {
    Row *x = &r->v[i];
    uint64_t k1 = ((uint64_t)(uint32_t)x->pi << 40)
                | ((uint64_t)(uint32_t)x->qs << 20) | (uint32_t)x->qe;
    if (!hset_add(&h1, k1 + 1)) continue;
    uint64_t k2 = ((uint64_t)(uint32_t)x->pi << 40)
                | ((uint64_t)(uint32_t)x->ts << 20) | (uint32_t)x->te;
    if (!hset_add(&h2, k2 + 1)) continue;
    r->v[keep++] = *x;
  }
  r->n = keep;
  free(h1.tab); free(h2.tab);
  return keep;
}

/* --------------------------------------------------------------------------
 * main
 * -------------------------------------------------------------------------- */
static void usage(void) {
  fprintf(stderr,
    "usage: mimic_scan --query Q.fa --db DB.fa --mode seeded|exhaustive --out O.tsv\n"
    "  [--kmer 4] [--seed-score-min 14] [--seed-gar-max 0.75] [--xdrop 6]\n"
    "  [--slack 20] [--gap-open 11] [--gap-ext 1] [--min-len 6] [--min-score 20]\n"
    "  [--n-iter 4] [--threads 4] [--dedup 1]\n");
  exit(2);
}

int main(int argc, char **argv) {
  const char *qpath = NULL, *dbpath = NULL, *modes = NULL, *outpath = NULL;
  Params P = {4, 14, 0.75, 6, 20, 11, 1, 6, 20, 4, 4};
  int do_dedup = 1;

  for (int i = 1; i < argc; i++) {
    const char *a = argv[i];
    if (!strcmp(a, "--query")) qpath = argv[++i];
    else if (!strcmp(a, "--db")) dbpath = argv[++i];
    else if (!strcmp(a, "--mode")) modes = argv[++i];
    else if (!strcmp(a, "--out")) outpath = argv[++i];
    else if (!strcmp(a, "--kmer")) P.kmer = atoi(argv[++i]);
    else if (!strcmp(a, "--seed-score-min")) P.seed_score_min = atoi(argv[++i]);
    else if (!strcmp(a, "--seed-gar-max")) P.seed_gar_max = atof(argv[++i]);
    else if (!strcmp(a, "--xdrop")) P.xdrop = atoi(argv[++i]);
    else if (!strcmp(a, "--slack")) P.slack = atoi(argv[++i]);
    else if (!strcmp(a, "--gap-open")) P.gap_open = atoi(argv[++i]);
    else if (!strcmp(a, "--gap-ext")) P.gap_ext = atoi(argv[++i]);
    else if (!strcmp(a, "--min-len")) P.min_len = atoi(argv[++i]);
    else if (!strcmp(a, "--min-score")) P.min_score = atoi(argv[++i]);
    else if (!strcmp(a, "--n-iter")) P.n_iter = atoi(argv[++i]);
    else if (!strcmp(a, "--threads")) P.threads = atoi(argv[++i]);
    else if (!strcmp(a, "--dedup")) do_dedup = atoi(argv[++i]);
    else usage();
  }
  if (!qpath || !dbpath || !modes || !outpath) usage();
  init_matrix();

  Recs q; memset(&q, 0, sizeof q);
  read_fasta(qpath, &q);
  if (q.n < 1) { fprintf(stderr, "empty query fasta\n"); return 1; }
  Recs prot; memset(&prot, 0, sizeof prot);
  read_fasta(dbpath, &prot);
  fprintf(stderr, "[io] query %s (%d aa), db %d sequences\n",
          q.v[0].acc, q.v[0].len, prot.n);

  Rows rows; memset(&rows, 0, sizeof rows);

  if (!strcmp(modes, "seeded")) {
    KIdx idx; memset(&idx, 0, sizeof idx);
    build_index(&idx, &prot, P.kmer);
    fprintf(stderr, "[io] k-mer index: %u entries\n", idx.nitems);
    scan_seeded(&prot, &q.v[0], &idx, &P, &rows);
    if (do_dedup) {
      int before = rows.n;
      int after = dedup(&rows);
      fprintf(stderr, "[dedup] %d -> %d rows\n", before, after);
    }
  } else if (!strcmp(modes, "exhaustive")) {
    int nt = P.threads;
    if (nt < 1) nt = 1;
    if (nt > 32) nt = 32;
    pthread_t th[32];
    ExhJob *jobs = (ExhJob *)calloc((size_t)nt, sizeof(ExhJob));
    for (int t = 0; t < nt; t++) {
      jobs[t].prot = &prot; jobs[t].query = &q.v[0]; jobs[t].P = &P;
      jobs[t].tid = t; jobs[t].nthreads = nt;
      memset(&jobs[t].rows, 0, sizeof(Rows));
    }
    for (int t = 1; t < nt; t++) pthread_create(&th[t], NULL, exh_worker, &jobs[t]);
    exh_worker(&jobs[0]);
    for (int t = 1; t < nt; t++) pthread_join(th[t], NULL);
    for (int t = 0; t < nt; t++)
      for (int i = 0; i < jobs[t].rows.n; i++) rows_push(&rows, &jobs[t].rows.v[i]);
  } else {
    fprintf(stderr, "unknown mode %s\n", modes);
    return 2;
  }

  FILE *o = fopen(outpath, "w");
  if (!o) { fprintf(stderr, "cannot write %s\n", outpath); return 1; }
  fprintf(o, "hit_id\tmode\tacc\tgene\tqstart\tqend\ttstart\ttend\taln_len\tn_id\t"
             "uscore\tgscore\tscore_core\tnid_core\tlen_core\tq_core_start\tq_core_end\t"
             "t_core_start\tt_core_end\textension_method\tgqseq\tgtseq\n");
  for (int i = 0; i < rows.n; i++) {
    Row *r = &rows.v[i];
    const Rec *sub = &prot.v[r->pi - 1];
    fprintf(o, "%s:%d:%d:%d:%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t",
            sub->acc, r->qs, r->qe, r->ts, r->te,
            r->mode ? "exhaustive" : "seeded",
            sub->acc, sub->gene, r->qs, r->qe, r->ts, r->te, r->alen, r->nid);
    if (r->uscore == NA_INT) fprintf(o, "NA\t"); else fprintf(o, "%d\t", r->uscore);
    if (r->gscore == NA_INT) fprintf(o, "NA\t"); else fprintf(o, "%d\t", r->gscore);
    fprintf(o, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n",
            r->cscore, r->cnid, r->clen, r->cqs, r->cqe, r->cts, r->cte,
            r->ext, r->gq, r->gt);
  }
  fclose(o);
  fprintf(stderr, "[out] %d rows -> %s\n", rows.n, outpath);
  return 0;
}
