/* ===========================================================================
 * pair_null.c — per-alignment composition-matched null model.
 *
 * Why this exists
 * ---------------
 * Every aggregate null model we ran (whole-proteome shuffles, windowed
 * shuffles, short-motif shuffles) answers a *library-level* question:
 * "does EBNA1 produce more / better hits than a composition-matched decoy?"
 * None of them answers the question a candidate list actually needs:
 *
 *      given the local composition of THIS EBNA1 segment and THIS human
 *      segment, how often would a random ordering of the same residues
 *      match this well?
 *
 * An aggregate p-value cannot be transferred to an individual pair.  This
 * tool computes the pair-level p directly, so every candidate in the final
 * list carries its own calibrated number instead of inheriting one from a
 * library-wide test (GPT review P0-2/P0-3; Biomni's E_emp has the same
 * defect, one level up).
 *
 * Null model
 * ----------
 * For a reported alignment (query span, target span):
 *   obs   = Smith-Waterman score of the real spans
 *   null  = SW score after independently shuffling BOTH spans, preserving
 *           each span's exact residue composition (Fisher-Yates).
 *   p     = (1 + #{null >= obs}) / (1 + B)
 *
 * Shuffling both sides destroys residue order on both sides while keeping
 * each side's composition, so the test asks the narrow question we can
 * actually answer: "is there anything here beyond composition?"  It is
 * deliberately the conservative direction -- a match that survives is not
 * an artefact of either partner's composition.
 *
 * Alignment kernel is byte-identical to mimic_scan.c (strict SW, BLOSUM62
 * 24x24, LINEAR gap = gap_open per residue, gap_ext never applied) so the
 * numbers are comparable with the reproduction runs.
 *
 * Build: gcc -O3 -march=native -o pair_null pair_null.c
 * Usage: pair_null --query Q.fa --db DB.fa --cands C.tsv --b 1000 \
 *                  --seed 42 --flank 5 --gap 11 --out OUT.tsv
 *
 * Cands TSV columns (header required, tab separated):
 *   target  qstart  qend  tstart  tend  alnlen  pident
 * (1-based inclusive spans, as emitted by mmseqs convertalis)
 * =========================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* --- BLOSUM62, identical to mimic_scan.c / R/blosum62.R ------------------- */
static const char AA24[] = "ARNDCQEGHILKMFPSTWYVBZX*";

static const int BL62[23][23] = {
{     4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0},
{    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1},
{    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1},
{    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1},
{     0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2},
{    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1},
{    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1},
{     0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1},
{    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1},
{    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1},
{    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1},
{    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1},
{    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1},
{    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1},
{    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2},
{     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0},
{     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0},
{    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2},
{    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1},
{     0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1},
{    -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1},
{    -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1},
{     0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1}
};

static int S[24][24];

static void init_matrix(void) {
  for (int i = 0; i < 23; i++)
    for (int j = 0; j < 23; j++) S[i][j] = BL62[i][j];
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

/* --- FASTA --------------------------------------------------------------- */
typedef struct { char acc[64]; char gene[64]; char *seq; unsigned char *enc; int len; } Rec;
typedef struct { Rec *v; int n, cap; } Recs;

static void recs_push(Recs *rs, Rec r) {
  if (rs->n == rs->cap) { rs->cap = rs->cap ? rs->cap * 2 : 256;
    rs->v = realloc(rs->v, (size_t)rs->cap * sizeof(Rec)); }
  rs->v[rs->n++] = r;
}

static void read_fasta(const char *path, Recs *rs) {
  FILE *f = fopen(path, "r");
  if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(2); }
  char *line = NULL; size_t cap = 0; ssize_t n;
  Rec cur; int have = 0; size_t seqcap = 0;
  while ((n = getline(&line, &cap, f)) > 0) {
    while (n > 0 && (line[n-1] == '\n' || line[n-1] == '\r')) line[--n] = 0;
    if (line[0] == '>') {
      if (have) recs_push(rs, cur);
      memset(&cur, 0, sizeof(cur));
      seqcap = 0; cur.seq = NULL; cur.enc = NULL;
      /* sp|ACC|NAME  or  ACC  ; gene from GN=... */
      const char *h = line + 1;
      const char *p1 = strchr(h, '|');
      if (p1) {
        const char *p2 = strchr(p1 + 1, '|');
        if (p2) { size_t L = (size_t)(p2 - p1 - 1); if (L > 63) L = 63;
                  memcpy(cur.acc, p1 + 1, L); cur.acc[L] = 0; }
      } else {
        size_t L = 0; while (h[L] && h[L] != ' ' && L < 63) L++;
        memcpy(cur.acc, h, L); cur.acc[L] = 0;
      }
      const char *gn = strstr(line, "GN=");
      if (gn) { const char *e = gn + 3; size_t L = 0;
                while (e[L] && e[L] != ' ' && L < 63) L++;
                memcpy(cur.gene, gn + 3, L); cur.gene[L] = 0; }
      have = 1;
    } else if (have) {
      for (ssize_t i = 0; i < n; i++) {
        char c = line[i];
        if (c == ' ' || c == '\t') continue;
        if ((size_t)cur.len + 1 >= seqcap) {
          seqcap = seqcap ? seqcap * 2 : 256;
          cur.seq = realloc(cur.seq, seqcap);
          cur.enc = realloc(cur.enc, seqcap);
        }
        cur.seq[cur.len] = c;
        cur.enc[cur.len] = (unsigned char)aaidx(c);
        cur.len++;
      }
    }
  }
  if (have) recs_push(rs, cur);
  free(line); fclose(f);
}

static Rec *find_rec(Recs *rs, const char *acc) {
  for (int i = 0; i < rs->n; i++)
    if (strcmp(rs->v[i].acc, acc) == 0) return &rs->v[i];
  return NULL;
}

/* --- strict Smith-Waterman, linear gap (identical to mimic_scan.c) -------- */
static int *g_Ha = NULL, *g_Hb = NULL, *g_tmp = NULL;
static int g_tcap = 0;

static int sw_score(const unsigned char *a, int la,
                    const unsigned char *b, int lb, int gap) {
  if (la <= 0 || lb <= 0) return 0;
  int need = lb + 2;
  if (need > g_tcap) {
    g_tcap = need;
    g_Ha = realloc(g_Ha, (size_t)g_tcap * sizeof(int));
    g_Hb = realloc(g_Hb, (size_t)g_tcap * sizeof(int));
    g_tmp = realloc(g_tmp, (size_t)g_tcap * sizeof(int));
  }
  int *Hp = g_Ha, *Hc = g_Hb;
  for (int j = 0; j < need; j++) { Hp[j] = 0; Hc[j] = 0; }
  int best = 0;
  for (int i = 1; i <= la; i++) {
    Hc[0] = 0;
    int ai = a[i-1];
    for (int j = 1; j <= lb; j++) {
      int d = Hp[j-1] + S[ai][b[j-1]];
      int u = Hp[j] - gap;
      int l = Hc[j-1] - gap;
      int v = d;
      if (u > v) v = u;
      if (l > v) v = l;
      if (v < 0) v = 0;
      Hc[j] = v;
      if (v > best) best = v;
    }
    int *t = Hp; Hp = Hc; Hc = t;
  }
  return best;
}

/* --- xorshift64* --------------------------------------------------------- */
static uint64_t g_state = 88172645463325252ULL;
static inline uint64_t xr(void) {
  g_state ^= g_state >> 12; g_state ^= g_state << 25; g_state ^= g_state >> 27;
  return g_state * 2685821657736338717ULL;
}
static void shuffle_inplace(unsigned char *v, int n) {
  for (int i = n - 1; i > 0; i--) {
    int j = (int)(xr() % (uint64_t)(i + 1));
    unsigned char t = v[i]; v[i] = v[j]; v[j] = t;
  }
}

/* --- sort helper for quantiles ------------------------------------------ */
static int cmp_int(const void *a, const void *b) {
  int x = *(const int *)a, y = *(const int *)b;
  return (x > y) - (x < y);
}

int main(int argc, char **argv) {
  const char *qpath = NULL, *dbpath = NULL, *cpath = NULL, *opath = NULL;
  int B = 1000, seed = 42, flank = 5, gap = 11;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--query") && i + 1 < argc) qpath = argv[++i];
    else if (!strcmp(argv[i], "--db") && i + 1 < argc) dbpath = argv[++i];
    else if (!strcmp(argv[i], "--cands") && i + 1 < argc) cpath = argv[++i];
    else if (!strcmp(argv[i], "--out") && i + 1 < argc) opath = argv[++i];
    else if (!strcmp(argv[i], "--b") && i + 1 < argc) B = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--seed") && i + 1 < argc) seed = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--flank") && i + 1 < argc) flank = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--gap") && i + 1 < argc) gap = atoi(argv[++i]);
    else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 2; }
  }
  if (!qpath || !dbpath || !cpath || !opath) {
    fprintf(stderr, "usage: pair_null --query Q.fa --db DB.fa --cands C.tsv --out O.tsv [--b 1000 --seed 42 --flank 5 --gap 11]\n");
    return 2;
  }
  if (B < 1) B = 1;
  g_state = 88172645463325252ULL ^ ((uint64_t)seed * 0x9E3779B97F4A7C15ULL);

  init_matrix();

  Recs q = {0}, db = {0};
  read_fasta(qpath, &q);
  if (q.n < 1) { fprintf(stderr, "query empty\n"); return 2; }
  fprintf(stderr, "[pair_null] query %s (%d aa)\n", q.v[0].acc, q.v[0].len);
  read_fasta(dbpath, &db);
  fprintf(stderr, "[pair_null] db %d proteins\n", db.n);

  FILE *cf = fopen(cpath, "r");
  if (!cf) { fprintf(stderr, "cannot open %s\n", cpath); return 2; }
  FILE *of = fopen(opath, "w");
  if (!of) { fprintf(stderr, "cannot write %s\n", opath); return 2; }

  fprintf(of, "target\tgene\tqstart\tqend\ttstart\ttend\talnlen\tpident"
              "\tobs_score\tnull_mean\tnull_sd\tnull_max\tnull_p95\tnull_p99"
              "\tn_ge\tb\tmargin_z\tp_pair\n");

  char *line = NULL; size_t lcap = 0; ssize_t n;
  int first = 1, n_done = 0;
  unsigned char *qs = NULL, *ts = NULL;
  int *nullv = malloc((size_t)B * sizeof(int));

  while ((n = getline(&line, &lcap, cf)) > 0) {
    if (first) { first = 0; if (strstr(line, "target") || strstr(line, "qstart")) continue; }
    while (n > 0 && (line[n-1] == '\n' || line[n-1] == '\r')) line[--n] = 0;
    if (n == 0) continue;
    char *fld[8]; int nf = 0;
    char *wp = line; fld[nf++] = wp;
    for (; *wp && nf < 8; wp++) if (*wp == '\t') { *wp = 0; fld[nf++] = wp + 1; }
    if (nf < 5) continue;
    const char *acc = fld[0];
    int qs1 = atoi(fld[1]), qe1 = atoi(fld[2]);
    int ts1 = atoi(fld[3]), te1 = atoi(fld[4]);

    Rec *t = find_rec(&db, acc);
    if (!t) { fprintf(stderr, "[warn] %s not in db\n", acc); continue; }

    int qlo = qs1 - 1, qhi = qe1;              /* 0-based half-open */
    int tlo = ts1 - 1 - flank, thi = te1 + flank;
    if (qlo < 0) qlo = 0;
    if (qhi > q.v[0].len) qhi = q.v[0].len;
    if (tlo < 0) tlo = 0;
    if (thi > t->len) thi = t->len;
    int lq = qhi - qlo, lt = thi - tlo;
    if (lq <= 0 || lt <= 0) continue;

    qs = realloc(qs, (size_t)lq);
    ts = realloc(ts, (size_t)lt);
    memcpy(qs, q.v[0].enc + qlo, (size_t)lq);
    memcpy(ts, t->enc + tlo, (size_t)lt);

    int obs = sw_score(qs, lq, ts, lt, gap);

    long sum = 0; double sumsq = 0.0; int mx = 0, nge = 0;
    for (int b = 0; b < B; b++) {
      shuffle_inplace(qs, lq);
      shuffle_inplace(ts, lt);
      int s = sw_score(qs, lq, ts, lt, gap);
      nullv[b] = s;
      sum += s;
      sumsq += (double)s * (double)s;
      if (s > mx) mx = s;
      if (s >= obs) nge++;
      /* restore the real residues for the next round */
      memcpy(qs, q.v[0].enc + qlo, (size_t)lq);
      memcpy(ts, t->enc + tlo, (size_t)lt);
    }
    qsort(nullv, (size_t)B, sizeof(int), cmp_int);
    int p95 = nullv[(int)(0.95 * (B - 1))];
    int p99 = nullv[(int)(0.99 * (B - 1))];
    double mean = (double)sum / (double)B;
    double var = (B > 1) ? (sumsq - (double)sum * (double)sum / (double)B) / (double)(B - 1) : 0.0;
    double sd = var > 0 ? sqrt(var) : 0.0;
    double z = sd > 0 ? (obs - mean) / sd : 0.0;
    double pval = (1.0 + (double)nge) / (1.0 + (double)B);

    fprintf(of, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.6g\n",
            acc, t->gene[0] ? t->gene : "-", qs1, qe1, ts1, te1,
            nf >= 6 ? atoi(fld[5]) : (qe1 - qs1 + 1),
            nf >= 7 ? fld[6] : "-",
            obs, mean, sd, mx, p95, p99, nge, B, z, pval);
    n_done++;
    if (n_done % 25 == 0) fprintf(stderr, "[pair_null] %d done\n", n_done);
  }

  free(line); free(nullv); free(qs); free(ts);
  fclose(cf); fclose(of);
  fprintf(stderr, "[pair_null] %d alignments scored, B=%d\n", n_done, B);
  return 0;
}
