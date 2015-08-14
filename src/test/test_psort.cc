#include <stdio.h>
#include <string.h>

#include "util.h"
#include "conc.h"

struct point {
  double x;
  double y;
};

struct point_less {
public:
  bool operator() (const point& a, const point& b) const {
    return a.x * a.x + a.y * a.y < b.x * b.x + b.y * b.y;
  }
};

int test_psort(bwt::idx_t n, bwt::idx_t m, uint64_t seed) {
  /* generate random points */
  unsigned short rg[3] = { 
    (unsigned short)((seed >> 32) & 65535), 
    (unsigned short)((seed >> 16) & 65535), 
    (unsigned short)((seed >>  0) & 65535), 
  };
  /* since this is a merge example, we must
     have at least 2 elements (1 element from
     each bwt to merge) */
  point * P = new point[n];

  for (bwt::idx_t j = 0; j < m; j++) {
    for (bwt::idx_t i = 0; i < n; i++) {
      P[i].x = erand48(rg);
      P[i].y = erand48(rg);
    }
    point_less lt;
    bwt::psort(P, P + n, lt, 10, 5);
    if (bwt::check_sorted(P, P + n, lt)) {
      printf("%3ld OK\n", j);
    } else {
      printf("%3ld NG\n", j);
      return 0;
    }
  }
  return 1;
}

int main(int argc, char ** argv) {
  size_t n = (argc > 1 ? atol(argv[1]) : 1000000);
  size_t m = (argc > 2 ? atol(argv[2]) : 100);
  uint64_t seed = (argc > 3 ? atol(argv[3]) : 20);

#if parallel_model == parallel_model_native_tbb
  int nw = 1;
  const char * nw_s = getenv("NW");
  if (nw_s) nw = atoi(nw_s);
  tbb::task_scheduler_init sched_init(nw);
#endif

  if (test_psort(n, m, seed)) {
    return 0;
  } else {
    return 1;
  }
}
