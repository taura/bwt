/* 
 * wavelet matrix
 */

#pragma once

#include <string.h>

#include "typedefs.h"
#include "util.h"
#include "conc.h"
#include "succinct_bit_vector.h"

/* 
   wavelet matrix for a string T efficiently
   supports rank(c, i) operation, returning
   the number of occurrence of c in L[0:i];
   in the literature it returns the number of
   occurrence up to and including L[i], but
   as a matter of following C convention
   we exclude L[i]

   usage:
   wavelet_matrix wm;
   wm.init(s, n, alpha0, alpha1);
   wm.rank(c, i);

   or more simply,

   wavelet_matrix wm(s, n, alpha0, alpha1);
   wm.rank(c, i);

   wm is a wavelet matrix for the string
   s[0:n], assuming all characters are in
   range [alpha0,alpha1).

   wm.rank(c, i) returns the number of times
   c occurs in s[0:i].

 */

namespace bwt {

  struct wavelet_matrix {
  private:
    alpha_t * s;
    idx_t n;	   /* number of characters */
    /* the range of characters (alpha0 ... alpha1,
       inclusive) */
    alpha_t alpha0;
    alpha_t alpha1;
    int alpha_width; /* number of bits per alphabets */
#if 0
    /* interval at which prefix sums are recorded */
    idx_t sum_interval;
    idx_t add_count_zero_gran;
#endif
    /* array of succinct bit vectors 
       bit_vectors[l] is the succinct bit vector for level l */
    succinct_bit_vector * bit_vectors;
    alpha_t * bitmaps;

  public:
    wavelet_matrix() {
      s = 0;
      bit_vectors = 0;
      bitmaps = 0;
    }

    /*  
	input t : array of n alphabets
	level l : the bit position with which to partition t
	output u : array to which the partitioned chars go
    */

  private:
    void partition_by_bit(alpha_t * t, alpha_t * u, int l,
			  idx_t sum_interval=10000,
			  idx_t add_zero_count_gran=10000) {
      int shift = alpha_width - l - 1;
      idx_t z = n;
      for (idx_t i = 0; i < n; i++) {
	z -= ((t[i] >> shift) & 1);
      }
      alpha_t mask0 = (1 << (alpha_width - l - 1)) - 1;
      alpha_t mask1 = (1 << alpha_width) - (1 << (alpha_width - l - 1));
      assert((mask0 & mask1) == 0);
      assert((mask0 | mask1) == (1 << alpha_width) - 1);

      /* we partition the n elements into segments,
	 each of which has sum_interval elements.
	 - we count 0/1 in each segment
	 - we calc the prefix sum of the two count arrays
      */
      idx_t ns = (n + sum_interval - 1) / sum_interval;
      idx_t (*offsets)[2] = (idx_t (*)[2])new_<idx_t>(ns * 2, "wavelet_matrix prefix sum array");

      /* count 0/1 in each segment of sum_interval characters */
      /* parallel (doall) */
      // for(idx_t si = 0l; si < ns; si++)
      pfor(idx_t, si, (idx_t)0, ns, (idx_t)1) {
	idx_t * counts = offsets[si];
	counts[0] = counts[1] = 0;
	idx_t ii = si * sum_interval; 
	for (idx_t i = ii; i < min(ii + sum_interval, n); i++) {
	  idx_t bit = ((t[i] >> shift) & 1);
	  counts[bit]++;
	}
      } end_pfor;
      /* take prefix sum */
      idx_t acc[2] = { 0, 0 };
      /* parallel (prefix sum) */
      for (idx_t si = 0; si < ns; si++) {
	idx_t * counts = offsets[si];
	idx_t t[2] = { counts[0], counts[1] };  /* tmp */
	/* set prefix sum */
	counts[0] = acc[0]; 
	counts[1] = acc[1];
	/* accumulate counts */
	acc[0] += t[0]; 
	acc[1] += t[1];
      }
      assert(0 <= acc[0]);
      assert(     acc[0] <= n);
      assert(0 <= acc[1]);
      assert(     acc[1] <= n);
      assert(acc[0] + acc[1] == n);
      /* add the number of zeros to all the ones' counts
	 to make them offset */
      /* parallel (doall) */
      // for(idx_t si = 0l; si < ns; si++)
      pfor(idx_t, si, (idx_t)0, ns, add_zero_count_gran) {
	offsets[si][1] += acc[0];
      } end_pfor;
      /* now offsets[x][0] is the index that
	 the first zero in the segment x 
	 should go and offsets[x][1] the 
	 index the first one in the segment x
	 should go */
      /* parallel (doall) */
      // for(idx_t si = 0l; si < ns; si++)
      pfor(idx_t, si, (idx_t)0, ns, (idx_t)1) {
	idx_t * offs = offsets[si];
	idx_t ii = si * sum_interval; 
	for (idx_t i = ii; i < min(ii + sum_interval, n); i++) {
	  alpha_t b = (t[i] >> shift) & 1;
	  assert(((b == 0) & (offs[b] < z)) |
		 ((b == 1) & (offs[b] < n)));
	  idx_t d = offs[b];
	  u[d] = (t[d] & mask1) | (t[i] & mask0);
	  offs[b] = d + 1;
	}
      } end_pfor;
      delete_((idx_t *)offsets, ns * 2, "wavelet_matrix prefix sum array");
    }

    /* input t : array of n characters; the l-th bit
       of them constitutes the level l wavelet matrix
       output u : address that can accommodate
       log(sigma) * n bits
    */
    void set_bit(alpha_t * u, idx_t bit_idx, alpha_t x) {
      idx_t i = bit_idx / 8;
      idx_t b = bit_idx % 8;
      volatile alpha_t * p = &u[i];
      for (long t = 0; t < 1000; t++) {
	alpha_t c = *p;
	/* set b-th bit of u[i] to x */
	if (cmp_and_set_8(p, c, c ^ ((-x ^ c) & (1 << b)))) {
	  return;
	}
      }
      fprintf(stderr, 
	      "failed to set bit %ld after 1000 retries. likely a bug\n",
	      bit_idx);
      check(0);
    }

    /* t is array of n chars, each char of
       which is alpha_width bits. transpose
       it into alpha_width bit sequences */

    void transpose_bits(alpha_t * t, alpha_t * u,
			idx_t gran=10000) {
      /* TODO: remove too many single byte writes */
      pfor(idx_t, i, (idx_t)0, n, gran) {
	alpha_t c = t[i];
	for (int l = 0; l < alpha_width; l++) {
	  /* set u's ((l * n) + i)-th bit */
	  int shift = alpha_width - l - 1;
	  set_bit(u, (l * n) + i, (c >> shift) & 1);
	}
      } end_pfor;
      pfor(int, l, 0, alpha_width, 1) {
	bit_vectors[l].init((uint8_t *)u, l * n, (l + 1) * n);
      } end_pfor;
    }

    /* the number of bits necessary to represent x
       x = 7 (0x111 ) -> 3 
       x = 8 (0x1000) -> 4
       it is l s.t. 2^(l-1) <= x < 2^l */
    int int_log2(alpha_t x_) {
      int l = 0;
      idx_t x = x_;
      idx_t y = 1;
      while (y <= x) {
	l++;
	y += y;
      }
      assert(x <  (1U << l));
      assert(x >= (1U << (l - 1)));
      return l;
    }

  public:
    /* initialize wavelet matrix for string starting
       from s and containing n alphabets (i.e., s[0:n]),
       assuming each character x is alpha0 <= x <= alpha1
       (note alpha1 is included) */
    void init(alpha_t * s_, idx_t n_, alpha_t * w,
	      alpha_t alpha0_=0, 
	      alpha_t alpha1_=255,
	      idx_t sum_interval=10000,
	      idx_t add_zero_count_gran=10000,
	      size_t memcpy_gran=10000,
	      idx_t transpose_gran=10000) {
      s = s_;
      n = n_;
      alpha0 = alpha0_;
      alpha1 = alpha1_;
      alpha_width = int_log2(alpha1);
      bit_vectors = new_<succinct_bit_vector>(alpha_width, "wavelet matrix bitvector array");
      bitmaps = new_<alpha_t>(n, "wavelet matrix bitmap array");
    
      assert(alpha_width > 0);
      assert(alpha_width <= 8);
      /* alternately use two arrays to sort them */
      alpha_t * t[2];
      t[ alpha_width      % 2] = w;
      t[(alpha_width + 1) % 2] = bitmaps;
      /* parallel (doall) */
      pmemcpy(t[0], s, sizeof(alpha_t) * n, memcpy_gran);

      for (int l = 0; l < alpha_width; l++) {
	/* move chars with l-th bit = 0 in front of the 
	   destination array. the first l bits in each 
	   char do not move */
	partition_by_bit(t[l % 2], t[(l + 1) % 2], l, 
			 sum_interval, add_zero_count_gran);
      }
      assert(bitmaps == t[(alpha_width + 1) % 2]);
      transpose_bits(t[alpha_width % 2], bitmaps, transpose_gran);
      assert(w == t[alpha_width % 2]);
    }

    void fini() {
      if (bit_vectors) {
	assert(bitmaps);
	for (int l = 0; l < alpha_width; l++) {
	  bit_vectors[l].fini();
	}
	delete_(bit_vectors, alpha_width, "wavelet matrix bitvector array");
	delete_(bitmaps, n, "workspace to build wavelet matrix 1");
      } else {
	assert(!bitmaps);
      }
    }

    /* return the number of times c appears in s[0:n] */

    __attribute__((optimize("unroll-loops")))
    idx_t rank(alpha_t c, idx_t i) {
      stat.start(ts_event_wm_rank);
      assert(0 <= i);
      assert(     i <= n);
      idx_t p = 0;
      idx_t q = i;
      for (int l = 0; l < alpha_width; l++) {

	/* invariant
	   bit_vectors[l].start(ts_event_[p,q) has the characters
	   whose first l bits match c and on the
	   left of i?
	*/
	int shift = alpha_width - l - 1;
	succinct_bit_vector& b = bit_vectors[l];
	if ((c >> shift) & 1) {
	  p = b.n_zeros + b.rank1(p);
	  q = b.n_zeros + b.rank1(q);
	} else {
	  p = b.rank0(p);
	  q = b.rank0(q);
	}
	assert(0 <= p);
	assert(     p <= n);
	assert(0 <= q);
	assert(     q <= n);
	assert(p <= q);
      }
      stat.end(ts_event_wm_rank);
      return q - p;
    }

    /* just for test */
    idx_t rank_slow(alpha_t c, idx_t i) {
      idx_t nc = 0;
      for (idx_t j = 0; j < i; j++) {
	if (s[j] == c) 
	  nc++;
      }
      return nc;
    }

  
  };

}
