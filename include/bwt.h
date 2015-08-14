#pragma once

#include <math.h>

#include "typedefs.h"
#include "opts.h"
#include "util.h"
#include "conc.h"
#include "suffix_array.h"
#include "wavelet_matrix.h"
#include "sampled_suffix_array.h"
#include "gap_array.h"

/* 

   how to build the bwt of a substring 
   and invert it back to the original
   substring.

   ...p xxxxxxxxe...
   ...x xxxxxxxe...
   ...x xxxxxxe...
   ...x xxxxxe...
   ...x xxxxe...
   ...x xxxe...
   ...x xxe...
   ...x xe...
   ...x e...
   
   sort suffixes 
   
      L
   ...x xxe...
   ...x xxxxxxe...
   ...x e...
   ...x xxxxe...
   ...x xxxxxe...
   ...p xxxxxxxxe...
   ...x xxxe...
   ...x xe...
   ...x xxxxxxxe...
   
   we obtained bwt in the L column.
   note that 
   (i) L contains a character p, which is not
   in the substring we consider 
   (ii) L does not contain a character e, 
   the last character of the substring.
   
   LF-mapping becomes more complex. 
   
        L SA
   0 ...x xxe...
   1 ...x xxxxxxe...
   2 ...x e...
   3 ...x xxxxe...
   4 ...x xxxxxe...
   5 ...p xxxxxxxxe...
   6 ...x xxxe...
   7 ...x xe...
   8 ...x xxxxxxxe...

   consider how to extend a row, say 

   6 ...x xxxe...

   by one character ahead to:

   3 ...x xxxxe...

   the usual LF-mapping is 

   LF(6) = the number of occurrences of chars < L[6] in the entire string
         + the number of occurrences of chars L[6] in L before position 6
   
   a natural extension of this is:
   
   LF(6) = the number of occurrences of chars < L[6] in the substring
         + the number of occurrences of chars L[6] in L before position 6

   but this is not correct, due to remarks
   (i) and (ii).

   recall that what we would like to get is
   the rank of suffix xxxxe in the suffix array.

   the first term just works, as long as we
   count the characters in the original
   substring, not L.  this can be done either
   by counting the original substring, or
   fixing L to account for the chracter p and
   e.

   the second term needs more careful
   treatment.

   this term 

   (a) counts the character p which
   should not be included;
   
   (b) and does not count the chracter e
   which should be included.
   
   in both cases, they are not issues if 
   p (or e) != L[6].
   
   (a) becomes an issue (only) if p = L[6].
   this is avoided by replacing the p with
   0 (the character we assume does not
   occur otherwise).
   
   (b) becomes an issue (only) if e = L[6].
   this can be fixed by adding one to the result
   obtained above if it is >= the position of e...
   in the suffix array, and L[6] = e.

   */

namespace bwt {

  struct bwt {

    /* given the original string T[0:n], this
       represents bwt of T[a:b] in L[a:b].
       since we L[a:b] is a permutation of
       T[a-1:b-1] and T[b-1] is not necessarily
       the end-of-string character (i.e. $),
       the structure also maintains auxiliary
       information as follows.

       (i) e_char: the character T[b-1] (not present in L[a:b])
       (ii) e_rank: the rank of T[b-1:]
    */
    const alpha_t * T;
    idx_t n;
    alpha_t * L;
    idx_t a;
    idx_t b;
    alpha_t e_char;		/* T[b-1] */

    /* sampled SA and related info */
    sampled_suffix_array ssa;
    idx_t s_rank;	    /* the rank of T[a:] */
    idx_t e_rank;	    /* the rank of T[b-1:] */

    idx_t * C;			/* C[c] = the number of chars < c in
				   the _original_ substring T[a:b] */
    wavelet_matrix wm;		/* wavelet matrix of L[a:b] - { T[a-1] } */

    bwt() { 
      C = 0; 
    }

    void init(const alpha_t * T_, idx_t n_, 
	      alpha_t * L_, idx_t a_, idx_t b_) {
      T = T_;
      n = n_;
      L = L_;
      a = a_;
      b = b_;
      e_char = T[b-1];

      s_rank = -1;
      e_rank = -1; 
      C = 0;
    }

    void init_extra(alpha_t * W, bwt_opt& opt) {
      count_alphabets(opt);
      build_wavelet_matrix(W, opt);
    }

    void fini(bwt_opt& opt) {
      stat.start(ts_event_bwt_fini);
      wm.fini();
      ssa.fini();
      if (C) delete_(C, opt.alpha_max + 1, "count");
      stat.end(ts_event_bwt_fini);
    }
  
    /* calculate the number of times characters
       smaller than each character occurs in
       the original substring T[a:b] (T[a],
       ..., T[b-1]).  since L[a:b] is a
       permutation of { T[a-1], T[a], ...  T[b
       - 2] }, we count L[a:b] and fix it, by
       decrementing the count of T[a-1] and 
       incrementing that of T[b-1].
       we do no longer have the original string
       T, but when constructing this bwt,
       T[a-1] must have been somewhere in L[a:b],
       and its position is recoreded in s_rank
       (that is, L[s_rank] = T[a-1]).  T[b-1] 
       is not in L[a:b], so it was explicitly
       saved in e_char. */

    void count_alphabets(bwt_opt& opt) {
      stat.start(ts_event_bwt_count_alphabets);
      assert(C == 0);
      /* memory allocation: can use a fixed address */
      C = new_<idx_t>(opt.alpha_max + 1, "count");
      for (idx_t c = 0; c <= opt.alpha_max; c++) {
	C[c] = 0;
      }
      /* parallel (doall with atomic increments) */
      //for(idx_t i = a; i < b; i++) 
      pfor(idx_t, i, a, b, opt.count_alpha_gran) {
	idx_t c = L[i];
	assert(c <= opt.alpha_max);
	atomic_inc(&C[c]);
      } end_pfor;
      assert(C[L[s_rank]] > 0);
      C[L[s_rank]]--;		/* decrement for T[a-1] */
      C[e_char]++;		/* increment for T[b-1] */
      /* prefix sum */
      idx_t t = 0;		/* total so far */
      for (idx_t c = 0; c <= opt.alpha_max; c++) {
	idx_t s = C[c];
	C[c] = t;
	t += s;
      }
      stat.end(ts_event_bwt_count_alphabets);
    }

    /* build wavelet matrix of L[a:b], excluding
       T[a-1] that should not be counted upon
       LF-mapping */
    void build_wavelet_matrix(alpha_t * W, bwt_opt& opt) {
      stat.start(ts_event_bwt_build_wavelet_matrix);
      /* assert wm has not been initialized before */
      alpha_t t = L[s_rank];
      /* temporarily 'hide' T[a-1] */
      L[s_rank] = 0;		
      /* build wavelet matrix */
      wm.init(L + a, b - a, W + a,
	      opt.alpha_min, 
	      opt.alpha_max, 
	      opt.wavelet_matrix_sum_interval, 
	      opt.wavelet_matrix_add_zero_count_gran,
	      opt.memcpy_gran,
	      opt.wavelet_matrix_transpose_gran);
      /* get T[a-1] back in place */
      L[s_rank] = t;
      stat.end(ts_event_bwt_build_wavelet_matrix);
    }

    /* see above for how this works  */
    idx_t LF_map_mod(alpha_t c, idx_t i) {
      stat.start(ts_event_bwt_LF_map_mod);
      idx_t k = C[c] + wm.rank(c, i - a) + a;
      k += ((k >= e_rank) & (c == e_char));
      stat.end(ts_event_bwt_LF_map_mod);
      return k;
    }

    /* see above for how this works  */
    idx_t LF_map_mod2(alpha_t c, idx_t i, idx_t j) {
      stat.start(ts_event_bwt_LF_map_mod2);
      idx_t k = C[c] + wm.rank(c, i - a) + a;
      if (c == e_char && k >= e_rank) {
	idx_less lt(T, n);
	if (lt(b - 1, j)) {
	  k++;
	}
      }
      stat.end(ts_event_bwt_LF_map_mod2);
      return k;
    }

    /* inverse bwt */
    void ibwt(alpha_t * I) {
      /* backward reconstruction */
      idx_t i = e_rank;
      alpha_t c = e_char;
      for (idx_t j = b; j > a; j--) {
	I[j - 1] = c;
	c = L[i];
	i = LF_map_mod(c, i);
      }
    }

    /* sample ns samples from sa[a:b] and put
       them into ssa.  we want to make sure a
       and b-1 are always sampled.  to this
       end, let d = (b - a - 1) / (ns - 1), and
       give high priorities to x's of which 
       (x - a) / d is close to an integer
    */
    int sample_sa(idx_t * sa, idx_t ns, bwt_opt& opt) {
    
      /* we put ns samples between
	 [a,b), (the range covered by t),
	 making sure a and b-1 are always
	 sampled.
       
	 a=a0 < a1 < ... < a_{ns-1} < b
       
	 j-th sample is placed at
       
	 i = a + (b - a) * j / ns
	 (here, / is an integer division)
       
	 that is, 
       
	 a + (b - a) * j / ns
	 <= i < a + (b - a) * j / ns + 1
       
	 (here, / is a real division)
       
	 (i - a - 1) * ns < (b - a) * j <= (i - a) * ns
       
      */
      stat.start(ts_event_bwt_sample_sa);
      assert(ns >= 2 || b - a == 1); /* for a and b - 1 */
      ssa.init(ns, opt.ssa_init_gran);
      /* parallel for (doall) */
      pfor(idx_t, r, a, b, opt.ssa_sample_gran) {
	/* T[i:]'s rank is r */
	idx_t i = sa[r];
	assert(a <= i);
	assert(     i <  b);
	if (i == a) {
	  assert(s_rank == -1);
	  s_rank = r;
	} 
	if (i + 1 == b) {
	  assert(e_rank == -1);
	  e_rank = r;
	}
	if (b - a == 1) {
	  ssa.put(r, i);
	} else {
	  /* TODO: avoid division */
	  idx_t x = floor((i - 1 - a) * (ns - 1) / (double)(b - a - 1));
	  idx_t y = floor((i     - a) * (ns - 1) / (double)(b - a - 1));
	  assert(x <= y);
	  if ((idx_t)x != (idx_t)y) {
	    assert((idx_t)x < (idx_t)y);
	    ssa.put(r, i);
	  } else {
	    assert(i != a     || b - a == 1);
	    assert(i != b - 1 || b - a == 1);
	  }
	}
      } end_pfor;
      assert(s_rank != -1);
      assert(e_rank != -1);
#if !defined(NDEBUG) || !NDEBUG
      if (opt.assert_level>=2) {
	assert(ssa.count() == min(ns, b - a));
      }
#endif
      ssa.n = min(ns, b - a);
      stat.end(ts_event_bwt_sample_sa);
      return 1;
    }

    /* take two sampled arrays and resample them so that 
       the merged arrays still have <= ns elements.
       we also guarantee a is always sampled */
    int resample_sa(sampled_suffix_array& l, sampled_suffix_array& r,
		    idx_t c, /* l.b = r.a */
		    idx_t ns, 
		    gap_array& gap,
		    bwt_opt& opt) {
      stat.start(ts_event_bwt_resample_sa);
      ssa.init(ns, opt.ssa_init_gran);

      /* we are going to choose ns samples,
	 out of ns0 samples. in the first i 
	 samples, we try to 


	 0 1  ...  nn-1  all samples

	 0    ...  ns-1  chosen samples

	 the i-th sample is chosen 
	 if (j <= i * (ns - 1) / (nn-1))
       

	 (i) pos=a and pos=b-1 must be chosen */
      idx_t nn = l.n + r.n;
      idx_t n_all_samples = 0;
      idx_t n_taken_samples = 0;

      /* calculate new ranks of left suffixes */
      /* paralle for */
      for (ssa_iterator it = l.pos_begin(); it.has_next(); it.next()) {
	/* s is an arbitrary sample in left */
	idx_t rank = it.e->rank, pos = it.e->pos;
	/* calc how many suffixes from right
	   are in front of me.  recall gap[i]
	   is the number of suffixes x from
	   right satisfying T[SA[i-1]:] < x <
	   T[SA[i]:], so, if I am the suffix
	   SA[rank], there are
	   (gap[a] + ... + gap[rank]) right suffix
	   in front of me.  since my original
	   rank in the left was rank, my rank
	   will become i + (gap[a] + ... + gap[rank]).
	*/
	assert(a <= rank);
	assert(     rank < c);
      
	if (n_taken_samples * (nn - 1) <= n_all_samples * (ns - 1)) {
	  idx_t g = gap.sum(rank);
	  assert(g <= b - c);
	  assert(a <= rank + g);
	  assert(     rank + g < b);
	  ssa.put(rank + g, pos);
	  if (pos == a) {
	    assert(s_rank == -1);
	    s_rank = rank + g;
	  }
	  n_taken_samples++;
	} else {
	  assert(pos != a);
	}
	n_all_samples++;
      }

      /* calc new rank of right suffixes */
      /* paralle for */
      for (ssa_iterator it = r.pos_begin(); it.has_next(); it.next()) {
	idx_t new_rank = it.e->rank, pos = it.e->pos;
	assert(a + c <= new_rank);
	assert(         new_rank < b + c);

	if (n_taken_samples * (nn - 1) <= n_all_samples * (ns - 1)) {
	  ssa.put(new_rank - c, pos);
	  if (pos == b - 1) {
	    assert(e_rank == -1);
	    e_rank = new_rank - c;
	  }
	  n_taken_samples++;
	} else {
	  assert(pos != b - 1);
	}
	n_all_samples++;
      }
      assert(s_rank != -1);
      assert(e_rank != -1);
      assert(ssa.n == ns || ssa.n == b - a);
      stat.end(ts_event_bwt_resample_sa);
      return 1;			/* OK */
    }

    /* compute SA[r] */
    idx_t sa(idx_t r) {
      stat.start(ts_event_bwt_sa);
      assert(ssa.n > 0);
      assert(a <= r);
      assert(     r < b);
      /* if this happens, it means the client gave
	 ns == 0 when building this. 
	 it is user's fault */
      idx_t x = 0;
      idx_t j = r;
      while (1) {
	assert(a <= j);
	assert(     j < b);
	assert(x < b - a);
	/* SA[j] + x = SA[i] */
	idx_t p = ssa.get(j);
	if (p != -1) {
	  assert(a <= p + x);
	  assert(     p + x < b);
	  stat.end(ts_event_bwt_sa);
	  return p + x;
	}
	j = LF_map_mod(L[j], j);
	x++;
      }
    }

    /* find i s.t. T[sa[i-1]:] < T[x:] < T[sa[i]:]
       two special cases
       when T[x:] < T[sa[a]:], return a
       when T[sa[b-1]:] < T[x:], return b

    */
    idx_t bin_search(idx_t x, bwt_opt& opt) {
      stat.start(ts_event_bwt_bin_search);
      idx_less lt(T, n);
      if (lt(x, sa(a))) {
	stat.end(ts_event_bwt_bin_search);
	return a;
      }
      if (lt(sa(b-1), x)) {
	stat.end(ts_event_bwt_bin_search);
	return b;
      }
      idx_t p = a, q = b - 1;
      /* T[p:] < s < T[q:] */
#if !defined(NDEBUG) || !NDEBUG
      if (opt.assert_level>=2) {
	assert(lt(sa(p), x));
	assert(lt(x, sa(q)));
      }
#else
      (void)opt;
#endif
      while (q - p > 1) {
#if !defined(NDEBUG) || !NDEBUG
	if (opt.assert_level>=2) {
	  assert(lt(sa(p), x));
	  assert(lt(x, sa(q)));
	}
#endif
	idx_t r = (p + q) / 2;
	assert(p < r);
	assert(r < q);
	if (lt(x, sa(r))) {
	  /* s < T[SA[r]:] */
	  q = r;
	} else {
	  p = r;
	}
      }
      /* T[SA[q-1]:] < T[x:] < T[SA[q]:] */
#if !defined(NDEBUG) || !NDEBUG
      if (opt.assert_level>=2) {
	assert(lt(sa(p), x));
	assert(lt(x, sa(q)));
      }
#endif
      assert(p == q - 1);
      stat.end(ts_event_bwt_bin_search);
      return q;
    }

  };

  /** input:
      \param T 
      \param n
      \param a
      \param b 
      \param SA
    
      T is an input string of n characters.
      T[a,b] is the substring we compute 
      the bwt for. (T[a:b] does include T[a] but
      not T[b]).
      SA[0:b-a] is the suffix array of T[a:b].
      it is a permutation of { a, a+1, ..., b-1 }.

      \param L 
      the result is written into L[a:b].
      L[a:b] becomes permutation of T[a-1:b-1].
      (T[-1] means T[n-1]).  
      note that L becomes a permutation of 
      { T[a-1], T[a], ..., T[b-2] }, missing
      T[b-1] originally present in T and introducing
      T[a-1] originally not present in T.
    
      this method also returns auxiliary
      information necessary to rebuild T[a:b]
      from L[a:b].  
    
      (i) T[b-1], the missing last character
      in the input substring.
      (ii) its poisition in the suffix array.
    
      note that in the usual bwt of an entire
      string T[b-1] is always '$', the end-of-string
      character and its position in the suffix
      array is always 0.
    
      since the reconstruction is done
      backward from the last character of the
      input string, we need to where the last
      suffix is. that is, the position x in
      the suffix array such that SA[x] = b -
      1.  note that in the regular bwt of the
      entire string, reconstruction

      \return the position in L at which the
      last suffix was stored (i.e., SA[x - a] ==
      b - 1)
  */
  bwt sa_to_bwt(const alpha_t * T, idx_t n, idx_t a, idx_t b, idx_t * SA, 
		alpha_t * L, alpha_t * W,
		bwt_opt& opt) {
    pfor(idx_t, r, a, b, opt.sa_to_bwt_gran) {
      idx_t i = SA[r];
      assert(a <= i);
      assert(     i <  b);
      /* remember where the smallest/largest suffixes went.
	 we also remember the 'lost character', the character
	 in the original substring but not included in L */
      if (i == 0) {
	L[r] = T[n - 1];
      } else {
	L[r] = T[i - 1];
      }
    } end_pfor;
    bwt bwt;
    bwt.init(T, n, L, a, b);
    /* sample suffix array */
    bwt.sample_sa(SA, min(bwt.b - bwt.a, opt.ssa_n_samples), opt);
    /* build wavelet matrix and char counts */
    bwt.init_extra(W, opt);
    return bwt;
  }

  /** compute bwt of T[a:b] and put the result
      into L[a:b].  return the position x in L
      to which the last suffix T[b-1:] was
      written. i.e., SA[x] = b - 1 */
  bwt bwt_leaf(const alpha_t * T, idx_t n, 
	       idx_t a, idx_t b, alpha_t * L, alpha_t * W,
	       bwt_opt& opt) {
    stat.start(ts_event_bwt_leaf);
    assert(a < b);
    /* memory allocation: can use a fixed address */
    idx_t * SA0 = new_<idx_t>(b - a, "leaf sa");
    idx_t * SA = SA0 - a;
    sa_range(T, n, a, b, SA, opt.sort_rec_threshold, opt.merge_rec_threshold);
    if (0) sa_show(T, n, a, b, SA);
    bwt res = sa_to_bwt(T, n, a, b, SA, L, W, opt);
    assert(SA[res.s_rank] == a);
    assert(SA[res.e_rank] + 1 == b);
    delete_(SA0, b - a, "leaf sa");
    stat.end(ts_event_bwt_leaf);
    return res;
  }

  /* how to merge two partial bwts (say, l
     and r). that is,
   
     l = bwt of T[a:b] (a permutation of T[a-1:b-1])
     r = bwt of T[b:c] (a permutation of T[b-1:c-1])
   
     we like to build 
   
     m = bwt of T[a:c] (a permutation of T[a-1:c-1])
   
     it amounts to finding a position to
     which each suffix of r should go in l.
     that is, for each i in [b,c), find j (in [a,b)) 
     s.t.
   
     T[SA[j]:] < T[i:] < T[SA[j+1]:]   (*)
   
     as special cases, 
   
     (i) if T[i:] is smaller than any suffix
     in l, we consider j to be a - 1.
   
     (ii) if T[i:] is larger than any suffix
     in l, we consider j to be b - 1.
   
     in other words, T[SA[a-1]:] is a hypothetical
     string smaller than any string and T[SA[b]:]
     a hypothetical larger than any string.
   
     the backbone of the algorithm is similar
     to reconstructing the original string
     from a bwt. we go backward from the
     last string (i.e., T[c-1:]).  the position
     of this substring is determined by any 
     search (e.g., binary search).
   
     once positions of suffixes up to a
     certain point, say, T[i:], we now want
     to determine the position of T[i-1:].
     this can be carried out in a manner
     simlar to LF-mapping.  Let's say the
     suffix array of l looks like below
   
     L SA
     0 ...x xxe...
     1 ...x xxxxxxe...
     2 ...x e...
     3 ...x xxxxe...
     4 ...x xxxxxe...
     5 ...p xxxxxxxxe...
     6 ...x xxxe...
     7 ...x xe...
     8 ...x xxxxxxxe...
   
     and T[i-1:] = cyyyyyyy...

     we assume that the poition of T[i:]
     (yyyyyyy...) is known. for example,
     it is known that T[3:] < T[i:] < T[4:].
   
     L SA
     0 ...x xxe...
     1 ...x xxxxxxe...
     2 ...x e...
     3 ...x xxxxe...
     c yyyyyyy...
     4 ...x xxxxxe...
     5 ...p xxxxxxxxe...
     6 ...x xxxe...
     7 ...x xe...
     8 ...x xxxxxxxe...
   
     now we ask where cyyyyyyy...
     lies between (i.e., find j such that
     T[SA[j]:] < T[i-1:] < T[SA[j+1]:];
     in the figure above, j = 3).
   
     this is exactly the modified LF-mapping.
  */

  /* calc gap array from l and r */

  /* build gap array when merging l and r.
     gap array is an array of which gap[i] 
     is the number of elements that fall
     between l[i] and l[i+1].

     we also construct a map describing
     which suffix in went to which gap.
     if suffix T[j] went to gap i, we
     record gap_map[j] = i
  */

  int build_gap(bwt& l, bwt& r, gap_array& gap, bwt_opt& opt) {
    stat.start(ts_event_bwt_build_gap);
    assert(l.T == r.T);
    assert(l.n == r.n);
    assert(l.L == r.L);
    assert(l.b == r.a);
    idx_less lt(l.T, l.n);
#if !defined(NDEBUG) || !NDEBUG
    if (opt.assert_level>=2) {
      assert(r.b - 1 == r.sa(r.e_rank));
    }
#endif
    r.ssa.sort_by_pos();
    /* parallel (build gap) */
    // for (idx_t begin = r.a; begin < r.b; begin += opt.build_gap_segment_sz)
    pfor_step(idx_t, begin, r.a, r.b, opt.build_gap_segment_sz) {
      idx_t g = -1;
      idx_t end = min(begin + opt.build_gap_segment_sz, r.b);
      ssa_reverse_iterator it = r.ssa.pos_reverse_begin_from(end - 1);
      for (idx_t j = end - 1; j >= begin; j--) {
	if (j == end - 1) {
	  /* find i s.t. T[SA[i-1]:] < T[b-1:] < T[SA[i]:] */
	  assert(g == -1);
	  g = l.bin_search(j, opt);
	} else {
	  /* given         T[SA[i-1]:]  < T[j:]   < T[SA[i]:] 
	     find i' s.t.  T[SA[i'-1]:] < T[j-1:] < T[SA[i']:] 
	     by LF-mapping */
	  assert(l.a <= g);
	  assert(       g <= l.b);
	  idx_t g_next = l.LF_map_mod2(r.T[j], g, j);
#if !defined(NDEBUG) || !NDEBUG
	  if (opt.assert_level>=2) {
	    assert(g_next == l.bin_search(j, opt));
	  }
#endif
	  g = g_next;
	}
	/* T[SA[g]:] < T[j:] < T[SA[g+1]:]  */
	assert(l.a <= g);
	assert(       g <= l.b);
#if !defined(NDEBUG) || !NDEBUG
	if (opt.assert_level>=2) {
	  assert(g == l.a || lt(l.sa(g - 1), j));
	  assert(g == l.b || lt(j, l.sa(g)));
	}
#endif
	gap.inc(g);			/* race condition */
	if (it.has_next() && it.c->pos == j) {
	  it.c->rank += g;
	  it.next();
	}
      }
    } end_pfor_step;
    gap.set_prefix_sum();
    stat.end(ts_event_bwt_build_gap);
    return 1;
  }

  /* given bwt of range T[a:b] and bwt of range T[b:c],
     merge them into a single bwt of range T[a:c] */
  bwt bwt_merge(bwt& l, bwt& r, alpha_t * w, bwt_opt& opt) {
    stat.start(ts_event_bwt_merge_bwt);
    assert(l.L == r.L);
    assert(l.b == r.a);
    gap_array gap;
    gap.init(l.a, l.b, opt.gap_sum_gran, opt.gap_init_gran);
    build_gap(l, r, gap, opt);
    /* gap[i] is the number of suffixes s from r s.t.
       T[SA[i-1]:] < s < T[SA[i]:]. therefore, the merged
       suffix array should look like:

       gap[a]   suffixes from r ; SA[a] ; 
       gap[a+1] suffixes from r ; SA[a+1] ;
       ...
       gap[b-1] suffixes from r ; SA[b-1] ;
       gap[b]   suffixes from r 
    */

    /* x/y/z : pointers to left/right/result array */
    stat.start(ts_event_bwt_merge_loop);
    /* parallel (doall) */
    // for(idx_t i = l.a; i < l.b + 1; l++) 
    pfor(idx_t, i, l.a, l.b + 1, opt.merge_gran) {
      /* gap[i] elements from right */
      idx_t e = r.a + gap.sum(i);
      for (idx_t j = e - gap.get(i); j < e; j++) {
	assert(j < r.b);
	assert(i - r.a + j < r.b);
	w[i - r.a + j] = r.L[j];
      }
      /* one element from left */
      if (i < l.b) {
	assert(i < l.b);
	assert(i - r.a + e < r.b);
	w[i - r.a + e] = l.L[i];
      }
    } end_pfor;
    stat.end(ts_event_bwt_merge_loop);
    stat.start(ts_event_bwt_memcpy);
    /* parallel (doall) */
    pmemcpy((void*)&l.L[l.a], (void*)&w[l.a], 
	    sizeof(alpha_t) * (r.b - l.a), opt.memcpy_gran);
    stat.end(ts_event_bwt_memcpy);
    bwt bwt;
    bwt.init(l.T, l.n, l.L, l.a, r.b);
    /* build samples for the merged array, based on samples
       from left and right */
    bwt.resample_sa(l.ssa, r.ssa, l.b, min(bwt.b - bwt.a, opt.ssa_n_samples), 
		    gap, opt);
    bwt.init_extra(w, opt);
    gap.fini();
    l.fini(opt);
    r.fini(opt);
    stat.end(ts_event_bwt_merge_bwt);
    return bwt;
  }

  /* recursively built bwt of T[a:b], which
     is a range in the entire string T of n bytes;
     the result goes to L[a:b]; W[a:b] is used
     as a scratch memory */

  bwt bwt_rec(const alpha_t * T, idx_t n, idx_t a, idx_t b, 
	      alpha_t * L, alpha_t * W,
	      bwt_opt& opt) {
    if (b - a <= opt.bwt_rec_threshold) {
      return bwt_leaf(T, n, a, b, L, W, opt);
    } else {
      idx_t c = (a + b) / 2;
      decl_task_group tg;
      bwt l;
      tg_run(tg, l = bwt_rec(T, n, a, c, L, W, opt));
      bwt r = bwt_rec(T, n, c, b, L, W, opt);
      tg_wait(tg);
      return bwt_merge(l, r, W, opt);
    }
  }

  bwt pmbwt(const alpha_t * T, idx_t n, alpha_t * L, bwt_opt& opt) {
    alpha_t * W = new_<alpha_t>(n, "workspace to merge");
    bwt t = bwt_rec(T, n, 0, n, L, W, opt);
    delete_(W, n, "workspace to merge");
    return t;
  }

  int check_equal(alpha_t * T, alpha_t * I, idx_t a, idx_t b) {
    for (idx_t i = a; i < b; i++) {
      if (!check(T[i] == I[i])) return 0;
    }
    return 1;			/* OK */
  }

  void random_init(alpha_t * L, idx_t n, unsigned short rg[3]) {
    for (idx_t i = 0; i < n; i++) {
      L[i] = (alpha_t)nrand48(rg);
    }
  }

}
