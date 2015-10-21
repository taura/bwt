#include <stdio.h>
#include <string.h>
#include <vector>

#include "util.h"
#include "bwt.h"

#define bwt_check(exp) bwt::check_(exp, #exp, __FILE__, __LINE__, __func__)

/* test concurrent put to ssa */
int test_ssa_conc(bwt::idx_t n, bwt::idx_t m, uint64_t seed) {
  /* generate random alphabets */
  unsigned short rg[3] = { 
    (unsigned short)((seed >> 32) & 65535), 
    (unsigned short)((seed >> 16) & 65535), 
    (unsigned short)((seed >>  0) & 65535), 
  };

  bwt::idx_t * V = new bwt::idx_t[n];
  bwt::bwt_opt opt;
  opt.set_defaults();
  bwt::mallocator mem(opt);

  for (bwt::idx_t j = 0; j < m; j++) {
    bwt::sampled_suffix_array ssa;
    ssa.init(n, mem);
    for (bwt::idx_t i = 0; i < n; i++) {
      V[i] = nrand48(rg);
    }
    // put n items, concurrently
#pragma omp parallel
    {
#pragma omp barrier      
#pragma omp for
      for (bwt::idx_t i = 0; i < n; i++) {
	ssa.put(V[i], V[i] + 1);
      }
    }
    // check
    for (bwt::idx_t i = 0; i < n; i++) {
      if (!bwt_check(ssa.get(V[i]) == V[i] + 1)) 
	return 0;
    }
    ssa.fini(mem);
    printf("%ld OK\n", j);
  }
  delete[] V;
  return 1;
}

int main(int argc, char ** argv) {
  size_t n = (argc > 1 ? atol(argv[1]) : 10000);
  size_t m = (argc > 2 ? atol(argv[2]) : 100);
  uint64_t seed = (argc > 3 ? atol(argv[3]) : 20);
  if (test_ssa_conc(n, m, seed)) {
    return 0;
  } else {
    return 1;
  }
}
