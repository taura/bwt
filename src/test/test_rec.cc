#include <stdio.h>
#include <string.h>

#include "util.h"
#include "conc.h"
#include "bwt.h"

#define bwt_check(exp) bwt::check_(exp, #exp, __FILE__, __LINE__, __func__)

int test_rec(bwt::idx_t n, bwt::idx_t m, 
	     bwt::alpha_t alpha0, bwt::alpha_t alpha1,
	     uint64_t seed) {
  /* generate random alphabets */
  assert(alpha0 > 0);
  unsigned short rg[3] = { 
    (unsigned short)((seed >> 32) & 65535), 
    (unsigned short)((seed >> 16) & 65535), 
    (unsigned short)((seed >>  0) & 65535), 
  };
  /* since this is a merge example, we must
     have at least 2 elements (1 element from
     each bwt to merge) */
  bwt_check(n >= 2);
  bwt::alpha_t * T = new bwt::alpha_t[n];
  for (bwt::idx_t i = 0; i + 1 < n; i++) {
    T[i] = alpha0 + (nrand48(rg) % (alpha1 - alpha0 + 1));
  }
  T[n - 1] = 0;
  
  /* bwt of T[a:b] into L[a:b] */
  bwt::alpha_t * L = new bwt::alpha_t[n];
  /* space to get the original string back from bwt */
  bwt::alpha_t * I = new bwt::alpha_t[n];

  for (bwt::idx_t j = 0; j < m; j++) {
    /* never forget to test the full string case */
    bwt::idx_t a = (j == 0 ? 0 : nrand48(rg) % (n - 1));
    assert(a < n - 1);
    bwt::idx_t b = (j == 0 ? n : a + 2 + nrand48(rg) % (n - 1 - a));
    assert(b <= n);
    assert(b - a >= 2);

    bwt::bwt_opt opt;
    opt.set_defaults(T, n);
    
    bwt::random_init(L, n, rg);
    /* bwt of the T[a:b] and check it's correct */
    printf("====== recursive ======\n");
    bwt::mstat.reset();
    bwt::mstat.level = 0;
    bwt::stat.reset();
    bwt::stat.level = 0;
    bwt::stat.start(bwt::ts_event_bwt_range_rec);
    bwt::tsc_t t0 = bwt::get_tsc();
    /* workspace for merge, necessary only for recursive algorithm */
    bwt::alpha_t * W = bwt::new_<bwt::alpha_t>(n, "workspace to merge");
    //bwt t = bwt_range_rec(T, n, a, b, threshold, ns, alpha0, alpha1, L, W);
    bwt::bwt t = bwt::bwt_rec(T, n, a, b, L, W, opt);
    bwt::delete_(W, n, "workspace to merge");
    bwt::tsc_t t1 = bwt::get_tsc();
    printf("%llu clocks to build bwt for %ld chars\n",
	   t1 - t0, b - a);

    bwt::stat.end(bwt::ts_event_bwt_range_rec);
    bwt::stat.print();
    bwt::random_init(I, n, rg);
    t.ibwt(I);
    t.fini(opt);
    bwt::check_equal(T, I, a, b);
    bwt::mstat.print();

    if (bwt::check_equal(T, I, a, b)) {
      printf("%lu OK\n", j);
    } else {
      printf("%lu NG\n", j);
      return 0;
    }
  }
  delete[] T;
  delete[] L;
  delete[] I;
  return 1;
}

int main(int argc, char ** argv) {
  size_t n = (argc > 1 ? atol(argv[1]) : 1000000);
  size_t m = (argc > 2 ? atol(argv[2]) : 100);
  bwt::alpha_t alpha0 = (argc > 4 ? (bwt::alpha_t)argv[4][0] : 1);
  bwt::alpha_t alpha1 = (argc > 5 ? (bwt::alpha_t)argv[5][0] : 255);
  uint64_t seed = (argc > 6 ? atol(argv[6]) : 20);

#if parallel_model == parallel_model_native_tbb
  int nw = 1;
  const char * nw_s = getenv("NW");
  if (nw_s) nw = atoi(nw_s);
  new tbb::task_scheduler_init(nw);
#endif

  if (test_rec(n, m, alpha0, alpha1, seed)) {
    return 0;
  } else {
    return 1;
  }
}
