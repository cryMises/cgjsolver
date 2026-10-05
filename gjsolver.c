#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef long long ll;

/* largest number of variables; above 10 the numbers overflow 64 bits for
 * ordinary inputs */
#define MAX 10
/* -1 sets its limits to the "bottom 2.5%" rows of its own steps and entries
 * tables (mu - 2 sigma), so they are known only after the seeds are run.
 * -4 has no seed statistics and keeps this fixed fit instead, fitted to the
 * average of the 10 rarest (lowest) results per random dense system,
 * coefficients -9..9, n = 1..10:
 *   steps   <= n^2 * (1 + n / STEPS_DIV)
 *   entries <= max(PEAK_FLOOR, 10^(PEAK_SLOPE * n - PEAK_SHIFT)) */
#define STEPS_DIV 20
#define PEAK_SLOPE 2.0
#define PEAK_SHIFT 6.5
#define PEAK_FLOOR 100
/* the two limits for the current n, set once n is read */
#define LIM_PEAK lim_peak
#define LIM_STEPS lim_steps
/* -1 always checks at least this many seeds */
#define MIN_TRIES 1000000u
/* -1 gives up after this many seeds */
#define MAX_TRIES 1000000u
/* -4 abandons a branch whose entries pass this; -2 and -3 use NEVER instead */
#define BIG 2000000000LL
/* no limit: only a real 64-bit overflow stops a branch */
#define NEVER LLONG_MAX
/* -2 shrinks a row only when an entry passes this, to avoid 64-bit overflow */
#define SHRINK_AT 100000000LL
/* search node limit for -2 */
#define NODE_LIMIT 8000000ULL
/* search node limit for -3 */
#define NODE_LIMIT_3 20000000ULL
/* search node limit for -4 */
#define NODE_LIMIT_4 10000000ULL
/* seconds of search for -2, -3 and -4 */
#define TIME_LIMIT 20
/* largest common denominator verify() will check before it skips */
#define VERIFY_CAP 1000000000000000000LL
/* most size-vs-steps plans kept by -4 */
#define FRONT_MAX 256
/* size of the step-count histogram in -1 */
#define HIST_MAX 8192
/* number of shrink rules tried by -4 */
#define NPOL 5
/* shrink thresholds of those rules */
static const ll POLICY[NPOL] = {0, 1000, 100000, 100000000LL,
                                4000000000000000000LL};
/* names of those rules, as printed in the -4 table */
static const char *PNAME[NPOL] = {"every row", "entries > 1e3", "entries > 1e5",
                                  "entries > 1e8", "never"};

/* rows of the steps and entries tables in the -1 table: mu-2sigma, mu-sigma,
 * mu, mu+sigma, mu+2sigma (labelled bottom 2.5% ... top 97.5%), where mu and
 * sigma are the mean and standard deviation of this run's own seeds (steps: of
 * the step counts; entries: of log10 of the largest entry, so the entry rows
 * are mu * 10^(k*sigma)) */
#define NFORM 5
#define NSTEPF 5
static const char *SIGL[NFORM] = {"bottom 2.5%", "bottom 16%", "middle 50%",
                                  "top 84%", "top 97.5%"};
/* width of the label column in the -1 tables */
#define LABW 28
/* -1 prints entries from this size up in scientific notation, smaller ones in
 * full */
#define SCI_FROM 1000000

static ll big_limit = BIG;

int n;
ll orig[MAX][MAX + 1];
ll a[MAX][MAX + 1];
int verbose = 1;
ll peak;
static double lim_peak;
static int lim_steps;

static void sci(char *b, size_t len, double v) {
  if (v < SCI_FROM)
    snprintf(b, len, "%.0f", v);
  else
    snprintf(b, len, "%.2e", v);
}

static double my_ln(double x) {
  double z = (x - 1) / (x + 1), z2 = z * z, t = z, sum = 0;
  for (int i = 1; i < 400; i += 2) {
    sum += t / i;
    t *= z2;
  }
  return 2 * sum;
}

static double my_exp(double t) {
  int k = 0;
  while (t > 0.5 || t < -0.5) {
    t /= 2;
    k++;
  }
  double term = 1, sum = 1;
  for (int i = 1; i < 20; i++) {
    term *= t / i;
    sum += term;
  }
  while (k--)
    sum *= sum;
  return sum;
}

static double my_pow(double x, double y) { return my_exp(y * my_ln(x)); }

/* 10^l as a long long, rounded and clamped so it cannot overflow */
static ll pow10_ll(double l) {
  double v = my_pow(10.0, l) + 0.5;
  return v >= 9.2e18 ? LLONG_MAX : (ll)v;
}

static double my_sqrt(double x) {
  if (x <= 0)
    return 0;
  double r = x > 1 ? x : 1;
  for (int i = 0; i < 100; i++)
    r = 0.5 * (r + x / r);
  return r;
}

/* log10 of a positive integer; the series in my_ln only converges quickly for
 * x < 10, so split off the decades first */
static double my_log10(ll v) {
  double m = (double)v;
  int e = 0;
  while (m >= 10) {
    m /= 10;
    e++;
  }
  return e + my_ln(m) / my_ln(10.0);
}

/* steps, largest entry and seed of every valid seed, so the -1 limits can be
 * set to the bottom-2.5% rows once the whole run is known */
static int st_all[MAX_TRIES];
static ll pk_all[MAX_TRIES];
static unsigned sd_all[MAX_TRIES];

static char *lim_peak_str(void) {
  static char b[24];
  sci(b, sizeof b, lim_peak);
  return b;
}

static char *obuf = NULL;
static size_t olen = 0, ocap = 0;

void say(const char *fmt, ...) {
  if (!verbose)
    return;
  va_list ap, ap2;
  va_start(ap, fmt);
  va_copy(ap2, ap);
  int need = vsnprintf(NULL, 0, fmt, ap);
  va_end(ap);
  if (need < 0) {
    va_end(ap2);
    return;
  }
  if (olen + (size_t)need + 1 > ocap) {
    size_t nc = ocap ? ocap : 4096;
    while (nc < olen + (size_t)need + 1)
      nc *= 2;
    obuf = realloc(obuf, nc);
    if (!obuf)
      exit(1);
    ocap = nc;
  }
  vsnprintf(obuf + olen, (size_t)need + 1, fmt, ap2);
  va_end(ap2);
  olen += (size_t)need;
}

void flush_out(void) {
  if (obuf) {
    fwrite(obuf, 1, olen, stdout);
    olen = 0;
  }
  fflush(stdout);
}

ll gcd(ll x, ll y) {
  x = llabs(x);
  y = llabs(y);
  while (y) {
    ll t = x % y;
    x = y;
    y = t;
  }
  return x;
}

void shuffle(int *v, int len) {
  for (int i = len - 1; i > 0; i--) {
    int j = rand() % (i + 1);
    int t = v[i];
    v[i] = v[j];
    v[j] = t;
  }
}

void print_matrix(void) {
  if (!verbose)
    return;
  int w = 4;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= n; j++) {
      char b[32];
      int l = snprintf(b, sizeof b, "%lld", a[i][j]);
      if (l > w)
        w = l;
    }
  for (int i = 0; i < n; i++) {
    say("   ");
    for (int j = 0; j < n; j++)
      say("%*lld ", w, a[i][j]);
    say("| %*lld\n", w, a[i][n]);
  }
}

void print_frac(ll num, ll den) {
  ll g = gcd(num, den);
  if (g == 0)
    g = 1;
  num /= g;
  den /= g;
  if (den < 0) {
    num = -num;
    den = -den;
  }
  if (den == 1)
    say("%lld", num);
  else
    say("%lld/%lld", num, den);
}

static int lin(ll m1, ll x, ll m2, ll y, ll *o) {
  ll t1, t2;
  if (__builtin_mul_overflow(m1, x, &t1) ||
      __builtin_mul_overflow(m2, y, &t2) || __builtin_sub_overflow(t1, t2, o))
    return 0;
  return 1;
}

static const char *verify(const ll *num, const ll *den) {
  __int128 L = 1;
  for (int j = 0; j < n; j++) {
    ll d = llabs(den[j]);
    ll g = gcd((ll)(L % d), d);
    L = L / g * d;
    if (L > (__int128)VERIFY_CAP)
      return "skipped";
  }
  for (int i = 0; i < n; i++) {
    __int128 lhs = 0;
    for (int j = 0; j < n; j++) {
      ll nj = den[j] < 0 ? -num[j] : num[j];
      ll dj = llabs(den[j]);
      lhs += (__int128)orig[i][j] * nj * (L / dj);
    }
    if (lhs != (__int128)orig[i][n] * L)
      return "FAILED";
  }
  return "ok";
}

void print_mult(ll num, ll den, int row) {
  ll g = gcd(num, den);
  if (g == 0)
    g = 1;
  num /= g;
  den /= g;
  if (den < 0) {
    num = -num;
    den = -den;
  }
  if (den == 1)
    say("%lldR%d", num, row);
  else
    say("(%lld/%lld)R%d", num, den, row);
}

typedef struct {
  ll m1, c2, mx;
} Op;

static int row_op(ll M[MAX][MAX + 1], int r, int pr, int c, Op *op) {
  ll p = M[pr][c], f = M[r][c], g = gcd(p, f);
  ll m1 = p / g, m2 = f / g;
  if (m1 < 0) {
    m1 = -m1;
    m2 = -m2;
  }
  ll mx = 0;
  for (int j = 0; j <= n; j++) {
    if (!lin(m1, M[r][j], m2, M[pr][j], &M[r][j]))
      return -2;
    if (llabs(M[r][j]) > mx)
      mx = llabs(M[r][j]);
  }
  op->m1 = m1;
  op->c2 = -m2;
  op->mx = mx;
  return 1;
}

static ll row_shrink(ll M[MAX][MAX + 1], int r) {
  ll rg = 0;
  for (int j = 0; j <= n; j++)
    rg = gcd(rg, M[r][j]);
  if (rg <= 1)
    return 1;
  for (int j = 0; j <= n; j++)
    M[r][j] /= rg;
  return rg;
}

static void print_op(int step, int r, int pr, const Op *op) {
  if (op->m1 == 1) {
    if (op->c2 > 0)
      say(" %d. R%d += (%lld)R%d\n", step, r + 1, op->c2, pr + 1);
    else
      say(" %d. R%d -= (%lld)R%d\n", step, r + 1, -op->c2, pr + 1);
  } else {
    if (op->c2 > 0)
      say(" %d. R%d = (%lld)R%d + (%lld)R%d\n", step, r + 1, op->m1, r + 1,
          op->c2, pr + 1);
    else
      say(" %d. R%d = (%lld)R%d - (%lld)R%d\n", step, r + 1, op->m1, r + 1,
          -op->c2, pr + 1);
  }
}

static int pivot_step(ll M[MAX][MAX + 1], int c, int pr, const int *rows, int m,
                      ll shrink_limit, int step0, ll *pk, int show) {
  int st = 0;
  for (int q = 0; q < m; q++) {
    int r = rows[q];
    if (M[r][c] == 0)
      continue;
    Op op;
    if (row_op(M, r, pr, c, &op) < 0)
      return -2;
    st++;
    if (pk && op.mx > *pk)
      *pk = op.mx;
    if (show) {
      print_op(step0 + st, r, pr, &op);
      print_matrix();
    }
    if (op.mx > shrink_limit) {
      ll rg = row_shrink(M, r);
      if (rg > 1) {
        st++;
        if (show) {
          say(" %d. ", step0 + st);
          print_mult(1, rg, r + 1);
          say("\n");
          print_matrix();
        }
      }
    }
  }
  return st;
}

static void init_run(void) {
  peak = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= n; j++) {
      a[i][j] = orig[i][j];
      if (llabs(a[i][j]) > peak)
        peak = llabs(a[i][j]);
    }
}

static int finish_run(const int *prow, int steps) {
  int pcol[MAX];
  for (int c = 0; c < n; c++)
    pcol[prow[c]] = c;
  for (int c = 0; c < n; c++) {
    int i = c;
    while (pcol[i] != c)
      i++;
    if (i != c) {
      for (int j = 0; j <= n; j++) {
        ll t = a[i][j];
        a[i][j] = a[c][j];
        a[c][j] = t;
      }
      int t = pcol[i];
      pcol[i] = pcol[c];
      pcol[c] = t;
      steps++;
      say(" %d. R%d <-> R%d\n", steps, c + 1, i + 1);
      print_matrix();
    }
  }
  say(" Diagonal\n");
  print_matrix();
  ll num[MAX], den[MAX];
  for (int i = 0; i < n; i++) {
    ll d = a[i][i];
    if (d != 1) {
      steps++;
      say(" %d. ", steps);
      print_mult(1, d, i + 1);
      say("\n");
    }
    num[i] = a[i][n];
    den[i] = d;
  }
  say("\n Solution\n");
  for (int i = 0; i < n; i++) {
    say("   x%d = ", i + 1);
    print_frac(num[i], den[i]);
    say("\n");
  }
  say("\n Steps: %d | Max entry: %lld | Check: %s\n", steps, peak,
      verify(num, den));
  return steps;
}

int solve(unsigned seed) {
  int steps = 0;
  srand(seed);
  init_run();
  say("Seed %u\n", seed);
  say(" Start\n");
  print_matrix();

  int order[MAX];
  for (int i = 0; i < n; i++)
    order[i] = i;
  shuffle(order, n);

  int used[MAX] = {0};
  int prow[MAX];

  for (int k = 0; k < n; k++) {
    int c = order[k];

    int cand[MAX], cnt = 0;
    for (int r = 0; r < n; r++)
      if (!used[r] && a[r][c] != 0)
        cand[cnt++] = r;
    if (cnt == 0)
      return -1;
    int pr = cand[rand() % cnt];
    used[pr] = 1;
    prow[c] = pr;
    say(" Pivot x%d (R%d)\n", c + 1, pr + 1);

    int rows[MAX], m = 0;
    for (int r = 0; r < n; r++)
      if (r != pr)
        rows[m++] = r;
    shuffle(rows, m);
    int st = pivot_step(a, c, pr, rows, m, 0, steps, &peak, 1);
    if (st < 0)
      return -2;
    steps += st;
  }
  return finish_run(prow, steps);
}

static const char *no_plan(void) {
  verbose = 0;
  int rc = solve(1);
  verbose = 1;
  if (rc == -1)
    return "Singular: no unique solution.";
  return big_limit == BIG ? "Numbers too large for the size cap (2e9)."
                          : "Numbers too large: 64-bit overflow.";
}

static int bestc[MAX], bestr[MAX], pathc[MAX], pathr[MAX];
static int bestsum;
static ll shrink_at = SHRINK_AT;

static ll apply_peak;
static ll bestpk;
static unsigned long long node_limit = NODE_LIMIT;
static clock_t deadline;

#define OUT_OF_TIME() (((nodes & 63) == 0) && clock() > deadline)
static unsigned long long nodes;
static int aborted;

static int apply_pivot(ll M[MAX][MAX + 1], int c, int r) {
  int rows[MAX], m = 0;
  for (int i = 0; i < n; i++)
    if (i != r)
      rows[m++] = i;
  int st = pivot_step(M, c, r, rows, m, shrink_at, 0, &apply_peak, 0);
  if (st < 0)
    return -1;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= n; j++)
      if (llabs(M[i][j]) > big_limit)
        return -1;
  return st;
}

typedef struct {
  int c, r, cost;
  ll ap;
} Cand;

static int cand_cmp(const void *x, const void *y) {
  const Cand *p = x, *q = y;
  if (p->cost != q->cost)
    return p->cost - q->cost;
  return (p->ap > q->ap) - (p->ap < q->ap);
}

static void dfs(int k, ll M[MAX][MAX + 1], int s, int *usedc, int *usedr) {
  if (aborted)
    return;
  if (k == n) {
    int rowof[MAX], total = s;
    for (int t = 0; t < n; t++)
      rowof[pathc[t]] = pathr[t];
    for (int c = 0; c < n; c++)
      if (M[rowof[c]][c] != 1)
        total++;
    int seen[MAX] = {0}, cycles = 0;
    for (int c = 0; c < n; c++) {
      if (seen[c])
        continue;
      cycles++;
      for (int x = c; !seen[x]; x = rowof[x])
        seen[x] = 1;
    }
    total += n - cycles;
    if (total < bestsum) {
      bestsum = total;
      for (int t = 0; t < n; t++) {
        bestc[t] = pathc[t];
        bestr[t] = pathr[t];
      }
    }
    return;
  }
  Cand cand[MAX * MAX];
  int nc = 0;
  for (int c = 0; c < n; c++) {
    if (usedc[c])
      continue;
    int nz = 0;
    for (int i = 0; i < n; i++)
      if (M[i][c] != 0)
        nz++;
    for (int r = 0; r < n; r++) {
      if (usedr[r] || M[r][c] == 0)
        continue;
      cand[nc].c = c;
      cand[nc].r = r;
      cand[nc].cost = nz - 1;
      cand[nc].ap = llabs(M[r][c]);
      nc++;
    }
  }
  qsort(cand, nc, sizeof(Cand), cand_cmp);
  for (int q = 0; q < nc; q++) {
    if (++nodes > node_limit || OUT_OF_TIME()) {
      aborted = 1;
      return;
    }
    ll N2[MAX][MAX + 1];
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= n; j++)
        N2[i][j] = M[i][j];
    int st = apply_pivot(N2, cand[q].c, cand[q].r);
    if (st < 0 || s + st >= bestsum)
      continue;
    pathc[k] = cand[q].c;
    pathr[k] = cand[q].r;
    usedc[cand[q].c] = 1;
    usedr[cand[q].r] = 1;
    dfs(k + 1, N2, s + st, usedc, usedr);
    usedc[cand[q].c] = 0;
    usedr[cand[q].r] = 0;
    if (aborted)
      return;
  }
}

int solve_plan(const int *pc, const int *pr_) {
  int steps = 0;
  init_run();
  say(" Start\n");
  print_matrix();
  int prow[MAX];
  for (int k = 0; k < n; k++) {
    int c = pc[k], pr = pr_[k];
    prow[c] = pr;
    say(" Pivot x%d (R%d)\n", c + 1, pr + 1);
    int rows[MAX], m = 0;
    for (int r = 0; r < n; r++)
      if (r != pr)
        rows[m++] = r;
    int st = pivot_step(a, c, pr, rows, m, shrink_at, steps, &peak, 1);
    if (st < 0)
      return -2;
    steps += st;
  }
  return finish_run(prow, steps);
}

static int cand_cmp3(const void *x, const void *y) {
  const Cand *p = x, *q = y;
  if (p->ap != q->ap)
    return (p->ap > q->ap) - (p->ap < q->ap);
  return p->cost - q->cost;
}

static void dfs3(int k, ll M[MAX][MAX + 1], int s, ll pk, int *usedc,
                 int *usedr) {
  if (aborted)
    return;
  if (k == n) {
    int rowof[MAX], total = s;
    for (int t = 0; t < n; t++)
      rowof[pathc[t]] = pathr[t];
    for (int c = 0; c < n; c++)
      if (M[rowof[c]][c] != 1)
        total++;
    int seen[MAX] = {0}, cycles = 0;
    for (int c = 0; c < n; c++) {
      if (seen[c])
        continue;
      cycles++;
      for (int x = c; !seen[x]; x = rowof[x])
        seen[x] = 1;
    }
    total += n - cycles;
    if (pk < bestpk || (pk == bestpk && total < bestsum)) {
      bestpk = pk;
      bestsum = total;
      for (int t = 0; t < n; t++) {
        bestc[t] = pathc[t];
        bestr[t] = pathr[t];
      }
    }
    return;
  }
  Cand cand[MAX * MAX];
  int nc = 0;
  for (int c = 0; c < n; c++) {
    if (usedc[c])
      continue;
    for (int r = 0; r < n; r++) {
      if (usedr[r] || M[r][c] == 0)
        continue;
      ll T[MAX][MAX + 1];
      for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
          T[i][j] = M[i][j];
      apply_peak = pk;
      int st = apply_pivot(T, c, r);
      if (st < 0)
        continue;
      cand[nc].c = c;
      cand[nc].r = r;
      cand[nc].cost = st;
      cand[nc].ap = apply_peak;
      nc++;
    }
  }
  qsort(cand, nc, sizeof(Cand), cand_cmp3);
  for (int q = 0; q < nc; q++) {
    if (cand[q].ap > bestpk)
      break;
    if (cand[q].ap == bestpk && s + cand[q].cost >= bestsum)
      continue;
    if (++nodes > node_limit || OUT_OF_TIME()) {
      aborted = 1;
      return;
    }
    ll N2[MAX][MAX + 1];
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= n; j++)
        N2[i][j] = M[i][j];
    apply_peak = pk;
    int st = apply_pivot(N2, cand[q].c, cand[q].r);
    pathc[k] = cand[q].c;
    pathr[k] = cand[q].r;
    usedc[cand[q].c] = 1;
    usedr[cand[q].r] = 1;
    dfs3(k + 1, N2, s + st, cand[q].ap, usedc, usedr);
    usedc[cand[q].c] = 0;
    usedr[cand[q].r] = 0;
    if (aborted)
      return;
  }
}

typedef struct {
  ll pk;
  int steps, pol;
  int c[MAX], r[MAX];
} Pt;
static Pt fr[FRONT_MAX];
static int nfr, curpol;

static int dominated(ll pk, int st) {
  for (int i = 0; i < nfr; i++)
    if (fr[i].pk <= pk && fr[i].steps <= st)
      return 1;
  return 0;
}

static void add_pt(ll pk, int st) {
  if (dominated(pk, st))
    return;
  int w = 0;
  for (int i = 0; i < nfr; i++)
    if (!(pk <= fr[i].pk && st <= fr[i].steps))
      fr[w++] = fr[i];
  nfr = w;
  if (nfr >= FRONT_MAX)
    return;
  fr[nfr].pk = pk;
  fr[nfr].steps = st;
  fr[nfr].pol = curpol;
  for (int t = 0; t < n; t++) {
    fr[nfr].c[t] = pathc[t];
    fr[nfr].r[t] = pathr[t];
  }
  nfr++;
}

static void dfs4(int k, ll M[MAX][MAX + 1], int s, ll pk, int *usedc,
                 int *usedr) {
  if (aborted || dominated(pk, s))
    return;
  if (k == n) {
    int rowof[MAX], total = s;
    for (int t = 0; t < n; t++)
      rowof[pathc[t]] = pathr[t];
    for (int c = 0; c < n; c++)
      if (M[rowof[c]][c] != 1)
        total++;
    int seen[MAX] = {0}, cycles = 0;
    for (int c = 0; c < n; c++) {
      if (seen[c])
        continue;
      cycles++;
      for (int x = c; !seen[x]; x = rowof[x])
        seen[x] = 1;
    }
    total += n - cycles;
    add_pt(pk, total);
    return;
  }
  Cand cand[MAX * MAX];
  int nc = 0;
  for (int c = 0; c < n; c++) {
    if (usedc[c])
      continue;
    for (int r = 0; r < n; r++) {
      if (usedr[r] || M[r][c] == 0)
        continue;
      ll T[MAX][MAX + 1];
      for (int i = 0; i < n; i++)
        for (int j = 0; j <= n; j++)
          T[i][j] = M[i][j];
      apply_peak = pk;
      int st = apply_pivot(T, c, r);
      if (st < 0)
        continue;
      cand[nc].c = c;
      cand[nc].r = r;
      cand[nc].cost = st;
      cand[nc].ap = apply_peak;
      nc++;
    }
  }
  qsort(cand, nc, sizeof(Cand), cand_cmp3);
  for (int q = 0; q < nc; q++) {
    if (dominated(cand[q].ap, s + cand[q].cost))
      continue;
    if (++nodes > node_limit || OUT_OF_TIME()) {
      aborted = 1;
      return;
    }
    ll N2[MAX][MAX + 1];
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= n; j++)
        N2[i][j] = M[i][j];
    apply_peak = pk;
    int st = apply_pivot(N2, cand[q].c, cand[q].r);
    pathc[k] = cand[q].c;
    pathr[k] = cand[q].r;
    usedc[cand[q].c] = 1;
    usedr[cand[q].r] = 1;
    dfs4(k + 1, N2, s + st, cand[q].ap, usedc, usedr);
    usedc[cand[q].c] = 0;
    usedr[cand[q].r] = 0;
    if (aborted)
      return;
  }
}

static int pt_cmp(const void *x, const void *y) {
  const Pt *p = x, *q = y;
  return (p->pk > q->pk) - (p->pk < q->pk);
}

static int trade_off(void) {
  ll pk0 = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= n; j++)
      if (llabs(orig[i][j]) > pk0)
        pk0 = llabs(orig[i][j]);
  nfr = 0;
  nodes = 0;
  node_limit = NODE_LIMIT_4;
  big_limit = BIG;
  int any_abort = 0;
  verbose = 0;
  for (int p = 0; p < NPOL; p++) {
    ll M[MAX][MAX + 1];
    int usedc[MAX] = {0}, usedr[MAX] = {0};
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= n; j++)
        M[i][j] = orig[i][j];
    shrink_at = POLICY[p];
    curpol = p;
    aborted = 0;
    deadline = clock() + (clock_t)(TIME_LIMIT / NPOL) * CLOCKS_PER_SEC;
    dfs4(0, M, 0, pk0, usedc, usedr);
    if (aborted)
      any_abort = 1;
  }
  verbose = 1;
  if (nfr == 0) {
    say("%s\n", any_abort ? "No plan found within the limit." : no_plan());
    flush_out();
    return 1;
  }
  qsort(fr, nfr, sizeof(Pt), pt_cmp);
  say("Size vs steps (%s)\n\n",
      any_abort ? "not proven, limit hit" : "proven best");
  say(" Max entry   Steps   Shrink\n");
  for (int i = 0; i < nfr; i++)
    say(" %9lld   %5d   %s\n", fr[i].pk, fr[i].steps, PNAME[fr[i].pol]);
  int pick = 0, found = 0;
  for (int i = 0; i < nfr; i++)
    if (fr[i].pk <= LIM_PEAK && (!found || fr[i].steps < fr[pick].steps)) {
      pick = i;
      found = 1;
    }
  if (found)
    say("\nShown: fewest steps with max entry <= %s\n\n", lim_peak_str());
  else
    say("\nShown: smallest numbers (none <= %s)\n\n", lim_peak_str());
  shrink_at = POLICY[fr[pick].pol];
  solve_plan(fr[pick].c, fr[pick].r);
  flush_out();
  return 0;
}

int plan_mode(int mode) {
  if (mode == 4)
    return trade_off();
  ll M[MAX][MAX + 1];
  int usedc[MAX] = {0}, usedr[MAX] = {0};
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= n; j++)
      M[i][j] = orig[i][j];
  bestsum = 1 << 30;
  nodes = 0;
  aborted = 0;
  verbose = 0;
  deadline = clock() + (clock_t)TIME_LIMIT * CLOCKS_PER_SEC;
  big_limit = NEVER;
  if (mode == 3) {
    ll pk0 = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= n; j++)
        if (llabs(M[i][j]) > pk0)
          pk0 = llabs(M[i][j]);
    shrink_at = 0;
    node_limit = NODE_LIMIT_3;
    bestpk = (ll)1 << 62;
    dfs3(0, M, 0, pk0, usedc, usedr);
  } else {
    shrink_at = SHRINK_AT;
    dfs(0, M, 0, usedc, usedr);
  }
  verbose = 1;
  if (bestsum == (1 << 30)) {
    say("%s\n", aborted ? "No plan found within the limit." : no_plan());
    flush_out();
    return 1;
  }
  if (mode == 3)
    say("Smallest numbers: max entry %lld, then %d steps\n%s (%llu nodes)\n\n",
        bestpk, bestsum, aborted ? "Not proven (limit hit)" : "Proven best",
        nodes);
  else
    say("Fewest steps: %d (n^2=%d)\n%s (%llu nodes)\n\n", bestsum, n * n,
        aborted ? "Not proven (limit hit)" : "Proven best", nodes);
  solve_plan(bestc, bestr);
  flush_out();
  return 0;
}

static void gap_steps(char *b, size_t len, int st) {
  if (st <= LIM_STEPS)
    snprintf(b, len, "within limit");
  else
    snprintf(b, len, "%d over", st - LIM_STEPS);
}

static void gap_peak(char *b, size_t len, ll pk) {
  if (pk <= LIM_PEAK)
    snprintf(b, len, "within limit");
  else
    snprintf(b, len, "%.3gx the limit", (double)pk / LIM_PEAK);
}

int main(void) {

  printf(
      " Input format:\n"
      "   n seed                 n = number of variables (1-%d)\n"
      "   a1 a2 ... an b         one line per equation: n coefficients, then "
      "the constant\n"
      " Integers only; use 0 for a missing variable. The system must have a"
      " unique solution.\n"
      " Example (2x+y=5, x-y=1):\n  2 1\n  2 1 5\n  1 -1 1\n"
      " Seed:\n"
      "  >0  replay that seed\n   0  random seed\n"
      "  -1  search seeds for the best run within the limits (bottom 2.5%% "
      "of steps and largest entries)\n"
      "  -2  fewest steps\n  -3  smallest numbers\n"
      "  -4  size vs steps table, then the best plan within the entries "
      "limit\n",
      MAX);
  fflush(stdout);

  ll seed_in;
  if (scanf("%d %lld", &n, &seed_in) != 2 || n < 1 || n > MAX || seed_in < -4 ||
      seed_in > 4294967295LL) {
    printf("Bad first line: need n seed (1<=n<=%d, seed -4..4294967295).\n",
           MAX);
    return 1;
  }
  lim_steps = (int)(n * n * (1.0 + (double)n / STEPS_DIV) + 0.5);
  lim_peak = my_pow(10.0, PEAK_SLOPE * n - PEAK_SHIFT);
  if (lim_peak < PEAK_FLOOR)
    lim_peak = PEAK_FLOOR;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= n; j++)
      if (scanf("%lld", &orig[i][j]) != 1) {
        printf("Bad input: integers only.\n");
        return 1;
      }

  {
    int ch;
    while ((ch = getchar()) != '\n' && ch != EOF) {
    }
  }

  if (seed_in <= -2)
    return plan_mode((int)-seed_in);

  unsigned seed;
  if (seed_in == -1) {
    int goal = n * n;
    static unsigned hist[HIST_MAX];
    int best = -1, worst = 0;
    double sum = 0;
    long ok = 0;
    int ok_steps = -1;
    ll ok_peak = 0;
    unsigned ok_seed = 0;
    unsigned best_seed = 1, tries = 0;
    unsigned goal_hit = 0, skipped = 0;
    unsigned nv = 0, steps_ok = 0, peak_ok = 0;
    double s_st = 0, s_st2 = 0, s_lp = 0, s_lp2 = 0;
    double psum = 0;
    ll pmin = 0, pmax = 0;
    ll fewest_peak = 0, lowest_peak = 0;
    unsigned lowest_seed = 0;
    int lowest_steps = 0;

    unsigned long long base = ((unsigned long long)time(NULL) * 2654435761ULL ^
                               (unsigned long long)clock() * 40503ULL) %
                              4294967295ULL;
    verbose = 0;
    for (unsigned s = 1; s <= MAX_TRIES; s++) {

      unsigned cur = (unsigned)((base + s) % 4294967295ULL) + 1;
      int st = solve(cur);
      tries++;
      if (st == -2) {
        skipped++;
        continue;
      }
      if (st < 0) {
        verbose = 1;
        say("Singular: no unique solution.\n");
        flush_out();
        return 1;
      }
      if (best < 0 || st < best) {
        best = st;
        best_seed = cur;
        fewest_peak = peak;
      }
      if (nv == 0 || peak < lowest_peak ||
          (peak == lowest_peak && st < lowest_steps)) {
        lowest_peak = peak;
        lowest_steps = st;
        lowest_seed = cur;
      }
      if (nv == 0 || peak < pmin)
        pmin = peak;
      if (peak > pmax)
        pmax = peak;
      psum += (double)peak;
      double lp = my_log10(peak > 0 ? peak : 1);
      st_all[nv] = st;
      pk_all[nv] = peak;
      sd_all[nv] = cur;
      s_st += st;
      s_st2 += (double)st * st;
      s_lp += lp;
      s_lp2 += lp * lp;
      nv++;
      if (st > worst)
        worst = st;
      sum += st;
      hist[st < HIST_MAX - 1 ? st : HIST_MAX - 1]++;
      if (st <= goal && !goal_hit)
        goal_hit = tries;

      if (tries >= MIN_TRIES && goal_hit)
        break;
    }
    verbose = 1;
    unsigned valid = tries - skipped;
    if (valid == 0) {
      say("Overflow: numbers too large for every seed.\n");
      flush_out();
      return 1;
    }
    char e1[24], e2[24], e3[24], g1[40], g2[40];
    double mu_st = s_st / valid, sd_st = my_sqrt(s_st2 / valid - mu_st * mu_st);
    double mu_lp = s_lp / valid, sd_lp = my_sqrt(s_lp2 / valid - mu_lp * mu_lp);
    const char *rule = "   ------------------------------------------------\n";

    /* limits = the bottom 2.5% rows of the two tables below */
    lim_steps = (int)(mu_st - 2 * sd_st + 1e-9);
    lim_peak = (double)pow10_ll(mu_lp - 2 * sd_lp);
    if (lim_peak < 1)
      lim_peak = 1;
    for (unsigned i = 0; i < nv; i++) {
      int sk = st_all[i] <= LIM_STEPS, pk = (double)pk_all[i] <= LIM_PEAK;
      steps_ok += sk;
      peak_ok += pk;
      if (sk && pk) {
        ok++;
        if (ok_steps < 0 || st_all[i] < ok_steps ||
            (st_all[i] == ok_steps && pk_all[i] < ok_peak)) {
          ok_steps = st_all[i];
          ok_peak = pk_all[i];
          ok_seed = sd_all[i];
        }
      }
    }

    /* none within both limits: take the seed closest to both, measuring each
     * overshoot in sigmas (steps by sigma of steps, entry by sigma of its
     * log10) and minimising the combined distance */
    if (ok == 0) {
      double lim_lp = mu_lp - 2 * sd_lp, dbest = 0;
      double ssd = sd_st > 0 ? sd_st : 1, lsd = sd_lp > 0 ? sd_lp : 1;
      for (unsigned i = 0; i < nv; i++) {
        double ds = (st_all[i] - lim_steps) / ssd;
        double dp = (my_log10(pk_all[i] > 0 ? pk_all[i] : 1) - lim_lp) / lsd;
        if (ds < 0)
          ds = 0;
        if (dp < 0)
          dp = 0;
        double d = ds * ds + dp * dp;
        if (ok_steps < 0 || d < dbest ||
            (d == dbest && (st_all[i] < ok_steps ||
                            (st_all[i] == ok_steps && pk_all[i] < ok_peak)))) {
          dbest = d;
          ok_steps = st_all[i];
          ok_peak = pk_all[i];
          ok_seed = sd_all[i];
        }
      }
    }

    say("\n SEARCH\n");
    say("   Seeds searched   %u", tries);
    if (skipped)
      say("  (%u skipped: overflow)", skipped);
    say("\n   Steps            min %d | mean %.1f | max %d\n", best,
        sum / valid, worst);
    sci(e1, sizeof e1, (double)pmin);
    sci(e2, sizeof e2, psum / valid);
    sci(e3, sizeof e3, (double)pmax);
    say("   Largest entry    min %s | mean %s | max %s\n", e1, e2, e3);

    say("\n STEPS PER SEED   (mean %.1f, sigma %.2f)\n", mu_st, sd_st);
    say("   %-12s %10s %10s %9s\n", "Group", "Steps <=", "Seeds", "% seeds");
    say("%s", rule);
    for (int k = 0; k < NSTEPF; k++) {
      int thr = (int)(mu_st + (k - 2) * sd_st + 1e-9);
      unsigned cum = 0;
      for (int v = 0; v < HIST_MAX && v <= thr; v++)
        cum += hist[v];
      say("   %-12s %10d %10u %8.2f%%\n", SIGL[k], thr, cum,
          100.0 * cum / valid);
    }

    say("\n LARGEST ENTRY PER SEED   (geometric mean 10^%.2f, sigma %.2f "
        "decades)\n",
        mu_lp, sd_lp);
    say("   %-12s %10s %10s %9s\n", "Group", "Entry <=", "Seeds", "% seeds");
    say("%s", rule);
    for (int k = 0; k < NFORM; k++) {
      double lthr = mu_lp + (k - 2) * sd_lp;
      ll pthr = pow10_ll(lthr);
      unsigned cum = 0;
      for (unsigned i = 0; i < nv; i++)
        if (pk_all[i] <= pthr)
          cum++;
      char t[24];
      sci(t, sizeof t, (double)pthr);
      say("   %-12s %10s %10u %8.2f%%\n", SIGL[k], t, cum, 100.0 * cum / valid);
    }

    say("\n LIMITS\n");
    say("   %-26s %10s %9s\n", "", "Seeds", "% seeds");
    say(" %s", rule + 1);
    char lab[48];
    snprintf(lab, sizeof lab, "steps <= %d", LIM_STEPS);
    say("   %-26s %10u %8.2f%%\n", lab, steps_ok, 100.0 * steps_ok / valid);
    snprintf(lab, sizeof lab, "entry <= %s", lim_peak_str());
    say("   %-26s %10u %8.2f%%\n", lab, peak_ok, 100.0 * peak_ok / valid);
    say("   %-26s %10ld %8.2f%%\n", "both", ok, 100.0 * ok / valid);

    say("\n SEEDS\n");
    gap_steps(g1, sizeof g1, best);
    gap_peak(g2, sizeof g2, fewest_peak);
    sci(e1, sizeof e1, (double)fewest_peak);
    say("   Fewest steps     seed %u\n", best_seed);
    say("                    %d steps (%s), max entry %s (%s)\n", best, g1, e1,
        g2);
    gap_steps(g1, sizeof g1, lowest_steps);
    gap_peak(g2, sizeof g2, lowest_peak);
    sci(e1, sizeof e1, (double)lowest_peak);
    say("   Smallest entry   seed %u\n", lowest_seed);
    say("                    %d steps (%s), max entry %s (%s)\n", lowest_steps,
        g1, e1, g2);
    if (ok > 0) {
      sci(e1, sizeof e1, (double)ok_peak);
      say("   Chosen           seed %u\n", ok_seed);
      say("                    %d steps, max entry %s (within limits)\n\n",
          ok_steps, e1);
      best_seed = ok_seed;
    } else {
      gap_steps(g1, sizeof g1, ok_steps);
      gap_peak(g2, sizeof g2, ok_peak);
      sci(e1, sizeof e1, (double)ok_peak);
      say("   Closest to both  seed %u\n", ok_seed);
      say("                    %d steps (%s), max entry %s (%s)\n\n", ok_steps,
          g1, e1, g2);
      best_seed = ok_seed;
    }
    seed = best_seed;
  } else {
    seed = (seed_in == 0) ? (unsigned)time(NULL) : (unsigned)seed_in;
  }

  int rc = solve(seed);
  if (rc == -2) {
    olen = 0;
    say("Overflow: numbers too large for this seed.\n");
    flush_out();
    return 1;
  }
  if (rc < 0) {
    olen = 0;
    say("Singular: no unique solution.\n");
    flush_out();
    return 1;
  }
  flush_out();
  return 0;
}
