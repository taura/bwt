/* 
 * psort.h
 *
 * SELIM G. AKL AND NICOLA SANTORO
 * Optimal Parallel Merging and Sorting Without Memory Conflicts 
 * IEEE transactions on computers C-36 11 1367-1369
 * 
 */
#pragma once

#include "typedefs.h"
#include "conc.h"

namespace bwt {

#if parallel_model == parallel_model_task

  template<typename T, typename Compare>
    idx_t choose_min_idx(T * l, idx_t begin, idx_t end, const Compare& lt) {
    idx_t m = begin;
    for (idx_t i = begin + 1; i < end; i++) {
      if (lt(l[i], l[m])) {
	m = i;
      }
    }
    return m;
  }

  /* generic version, which workks whether d and l point to the same vector */
  template<typename T, typename Compare>
    inline void ins_sort_range(T * l, T * r, T * d, const Compare& lt) {
    idx_t n = r - l;
    for (idx_t i = 0; i < n; i++) {
      idx_t j = choose_min_idx(l, i, n, lt);
      T t = l[i];
      d[i] = l[j];
      l[j] = t;
    }
  }

  /* 
     binary search piv in [a,b). more precisely,
     find p s.t.
     (1) p == b || piv <= p[0], and
     (2) p == a || p[0] < piv
  */
  template<typename T, typename Compare>
    T * find_idx(T * a_beg, T * a_end, T piv, const Compare& lt) {
    if (a_beg == a_end)       return a_beg;
    if (!(lt(a_beg[0], piv))) return a_beg;
    if (lt(a_end[-1], piv))   return a_end;
    T * p = a_beg;
    T * q = a_end - 1;
    assert(lt(p[0], piv));
    assert(!(lt(q[0], piv)));
    while (q - p > 1) {
      assert(lt(p[0], piv));
      assert(!(lt(q[0], piv)));
      T * r = p + (q - p) / 2;
      if (lt(r[0], piv)) {
	p = r;
      } else {
	q = r;
      }
    }
    return q;
  }

  /* merge [a,b) and [c,d) into [t,..) */
  template<typename T, typename Compare>
    void p_merge(T * a_beg, T * a_end, 
		 T * b_beg, T * b_end, 
		 T * t_beg, const Compare& lt, 
		 idx_t merge_rec_threshold) {
    if ((a_end - a_beg) + (b_end - b_beg) < merge_rec_threshold) {
      /* serial merge */
      T * p = a_beg;
      T * q = b_beg;
      T * r = t_beg;
      while (p < a_end && q < b_end) {
	if (lt(*p, *q)) *r++ = *p++;
	else *r++ = *q++;
      }
      while (p < a_end) *r++ = *p++;
      while (q < b_end) *r++ = *q++;
    } else {
      T * p, * q;
      if (a_end - a_beg > b_end - b_beg) {
	p = a_beg + (a_end - a_beg) / 2;
	q = find_idx(b_beg, b_end, *p, lt);
      } else {
	q = b_beg + (b_end - b_beg) / 2;
	p = find_idx(a_beg, a_end, *q, lt);
      }
      idx_t n1 = p - a_beg;
      idx_t n2 = q - b_beg;
      decl_task_group tg;
      /* merge [a_beg,p) and [b_beg,q) into t */
      /* merge [p,b) and [q,d) into t + (p - a) + (q - c) */
      tg_run(tg, { 
	  p_merge(a_beg, p, b_beg, q, t_beg,           lt, merge_rec_threshold); 
	});
      p_merge(p, a_end, q, b_end, t_beg + n1 + n2, lt, merge_rec_threshold);
      tg_wait(tg);
    }
  }

  /* cilk sort [a_beg,a_end) into s, using t as temporary */
  template<typename T, typename Compare>
    void sort_range(T * a_beg, T * a_end, T * t_beg, const Compare& lt, 
		    int dest, idx_t sort_rec_threshold, idx_t merge_rec_threshold) {
    idx_t n = a_end - a_beg;
    if (n <= sort_rec_threshold) {
      ins_sort_range(a_beg, a_end, (dest == 0 ? a_beg : t_beg), lt);
    } else {
      idx_t nh = n / 2;
      T * c = a_beg + nh;
      decl_task_group tg;
      tg_run(tg, { 
	  sort_range(a_beg, c, t_beg, lt,
		     1 - dest, sort_rec_threshold, merge_rec_threshold);
	});
      sort_range(c, a_end, t_beg + nh, lt,
		 1 - dest, sort_rec_threshold, merge_rec_threshold);
      tg_wait(tg);
      T * s = (dest == 0 ? t_beg : a_beg);
      T * d = (dest == 0 ? a_beg : t_beg);
      p_merge(s, s + nh, s + nh, s + n, d, lt, merge_rec_threshold);
    }
  }

  template<typename T, typename Compare>
    void task_parallel_sort(T * a_beg, T * a_end, const Compare& lt, 
			    mallocator& mem, mem_reason_kind_t reason,
			    idx_t sort_rec_threshold=30, 
			    idx_t merge_rec_threshold=1000) {
    idx_t n = a_end - a_beg;
    T * t_beg = mem.new_<T>(n, reason);
    sort_range(a_beg, a_end, t_beg, lt, 0,
	       sort_rec_threshold, merge_rec_threshold);
    mem.delete_(t_beg, n, reason);
  }

#endif

  template<typename T, typename Compare>
    int check_sorted(T * a, T * e, const Compare& lt) {
    idx_t n = e - a;
    for (idx_t i = 0; i < n - 1; i++) {
      if (!check(lt(a[i], a[i + 1]))) return 0; /* NG */
    }
    return 1;			/* OK */
  }

}
