/* 
 * data.h
 */
#pragma once
#include <stdio.h>

alpha_t * gen_random_string(idx_t n, alpha_t alpha0, alpha_t alpha1, 
			    unsigned short rg[3]) {
  alpha_t * T = new alpha_t[n];
  for (idx_t i = 0; i + 1 < n; i++) {
    T[i] = alpha0 + (nrand48(rg) % (alpha1 - alpha0 + 1));
  }
  T[n - 1] = 0;
  return T;
}

alpha_t * gen_random_string_s(idx_t n, alpha_t alpha0, alpha_t alpha1, 
			      unsigned long seed) {
  unsigned short rg[3] = { (seed >> 32) & ((1 << 16) - 1), 
			   (seed >> 16) & ((1 << 16) - 1), 
			   (seed >>  0) & ((1 << 16) - 1) };
  return gen_random_string(n, alpha0, alpha1, rg);
}

alpha_t * read_string_from_file(idx_t n, const char * filename) {
  FILE * fp = fopen(filename, "rb");
  if (!fp) die("fopen");
  if (n == 0) {
    if (fseek(fp, 0, SEEK_END) != 0) die("fseek");
    n = ftell(fp);
    rewind(fp);
  }
  alpha_t * T = new alpha_t[n];
  idx_t m = fread(T, sizeof(alpha_t), n, fp);
  assert(m == n);
  T[n] = 0;
  return T;
}
