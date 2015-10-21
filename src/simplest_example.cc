#include "bwt.h"

int main() {
  /* input string (must end with null character) */
  const bwt::alpha_t * T = (bwt::alpha_t *)"mississippi";
  /* length must include the terminating null character */
  int n = strlen((const char *)T) + 1;
  /* allocate space for BWT */
  bwt::alpha_t * L = new bwt::alpha_t[n];

  /* set various parameters (tuning knobs) */
  bwt::bwt_opt opt;
  opt.set_defaults(T, n);
  bwt::mallocator mem(n, opt);

  /* do the real work */
  bwt::pmbwt(T, n, L, 0, mem, opt);

  printf("input: %s\n", T);
  printf("bwt:   ");
  for (int i = 0; i < n; i++) putchar(L[i]);
  putchar('\n');
  return 0;
}
/*
  mississippi$
  ississippi$m
  ssissippi$mi
  sissippi$mis
  issippi$miss
  ssippi$missi
  sippi$missis
  ippi$mississ
  ppi$mississi
  pi$mississip
  i$mississipp
  $mississippi

sort:

  $mississippi
  i$mississipp
  ippi$mississ
  issippi$miss
  ississippi$m
  mississippi$
  pi$mississip
  ppi$mississi
  sippi$missis
  sissippi$mis
  ssippi$missi
  ssissippi$mi

 */
