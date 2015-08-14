#pragma once

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "typedefs.h"
#include "util.h"
#include "conc.h"

namespace bwt {

  /* given i and j, lexicographically compare T[i] and T[j] */
  struct idx_less {
  idx_less(const alpha_t * T, idx_t n) 
  : T(T), n(n) {}
    const alpha_t * T;
    idx_t n;
  public:
    bool operator() (const idx_t& a, const idx_t& b) const {
      for (idx_t i = a, j = b; ; i++, j++) {
	assert(i < n);
	assert(j < n);
	if (T[i] < T[j]) return 1;
	if (T[i] > T[j]) return 0;
      }
    }
  };

  /* input:
     T[a:b] input string
     output:
     SA[0:b-a] suffix array of T[a:b] */
  int sa_range(const alpha_t * T, idx_t n, 
	       idx_t a, idx_t b, idx_t * SA,
	       idx_t sort_rec_threshold=30, 
	       idx_t merge_rec_threshold=1000) {
    idx_less lt(T, n);
    for (idx_t i = a; i < b; i++) {
      SA[i] = i;
    }
    /* parallel sort? */
    psort(SA + a, SA + b, lt, sort_rec_threshold, merge_rec_threshold);
    return 1;
  }

  /* given T and its suffix array for T[a:b], 
     print suffixes in the lexicographical order.
     each line is a cyclic shift, the suffix
     followed by the original string up to 
     one character before the suffix. so you
     will sse something like:
   
     L F    L
     o 
     h ello
     $ hello
     e llo
     l lo
     l o
  */

  void sa_show(const alpha_t * T, idx_t n, idx_t a, idx_t b, idx_t * SA) {
    alpha_t * C = new alpha_t[n];
    memcpy(C, T, sizeof(alpha_t) * n);
    for (idx_t i = a; i < b; i++) {
      assert(i >= a);
      assert(i < b);
      idx_t o = SA[i];
      assert(o >= a);
      assert(o < b);
      alpha_t p = (o ? C[o - 1] : C[n - 1]);
      printf("%2lu %2lu %c %s\n", i, o, (p ? p : '$'), C + o);
    }
    delete C;
  }


}

