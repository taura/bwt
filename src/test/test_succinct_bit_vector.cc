#include <stdio.h>
#include "succinct_bit_vector.h"

int test_succinct_bit_vector(bwt::idx_t n, bwt::idx_t m) {
  unsigned short rg[3] = { 1, 2, 3 };
  bwt::idx_t begin = 6;		/* begin bit position */
  bwt::idx_t end = begin + n;	/* end bit position */
  /* n bits -> nbytes */
  bwt::idx_t nb = (end + 7) / 8;	/* number of bytes */
  bwt::idx_t offset = 3;
  uint8_t * a0 = new uint8_t[nb + offset];
  uint8_t * a = a0 + offset;
  for (bwt::idx_t i = 0; i < nb; i++) {
    a[i] = (uint8_t)(nrand48(rg) % 256);
  }
  bwt::succinct_bit_vector s;
  bwt::bwt_opt opt;
  opt.set_defaults(0, n);
  bwt::mallocator mem(n, opt);
  
  printf("building a succinct bit vector bits@%p[%ld,%ld] (%ld bits) ...\n", 
	 a, begin, end, (end - begin));
  bwt::tsc_t c0 = bwt::get_tsc();
  s.init(a, begin, end, mem);
  bwt::tsc_t c1 = bwt::get_tsc();
  printf("took %llu cycles\n", c1 - c0);
  for (bwt::idx_t i = 0; i < m; i++) {
    bwt::idx_t p = (i == 0 ? end - begin : (bwt::idx_t)nrand48(rg) % (end - begin + 1));
    bwt::tsc_t c0 = bwt::get_tsc();
    bwt::idx_t x = s.rank1(p);
    bwt::tsc_t c1 = bwt::get_tsc();
    bwt::idx_t y = s.rank_slow(0, p);
    if (x == y) {
      printf("[OK] %llu cycles to find rank(%lu) = %lu\n", 
	     c1 - c0, p, x);
    } else {
      printf("[NG] %llu cycles to find rank(%lu) = %lu != %lu\n", 
	     c1 - c0, p, x, y);
      return 0;
    }
  }
  s.fini(mem);
  delete a0;
  printf("OK\n");
  return 1;			/* OK */
}

int main(int argc, char ** argv) {
  size_t n = (argc > 1 ? atol(argv[1]) : 1000000);
  size_t m = (argc > 2 ? atol(argv[2]) : 100);
  if (test_succinct_bit_vector(n, m)) {
    return 0;
  } else {
    return 1;
  }
}
