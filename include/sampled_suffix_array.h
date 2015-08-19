/* 
 * ssa.h
 */
#pragma once

#include "typedefs.h"

/* sampled suffix array is conceptually
   a sparse array of (rank, pos).
   in order to facilitate merging two
   sampled arrays, however, it supports
   a special iterator, called reverse_pos_iterator,
   to iterate all entries in a 

   sampled_sa ssa;

   (i) ssa.put(rank, pos);
   (ii) ssa.get(rank);
   (iii) usual iterator
   for (it = ssa.begin(); it != ssa.end(); it++) {
     ...
   }
   (iv) position-based operator
   for (it = ssa.reverse_pos_begin();
        it != ssa.reverse_pos_end();
        it++) {
     it->rank ...;
     it->pos ...;
     it->gap = x;
   }
   
 */

namespace bwt {

  struct ssa_entry {
    idx_t rank;			/* key */
    idx_t pos;			/* value */
  };

  struct ssa_pos_lt {
    int operator() (const ssa_entry& a, const ssa_entry& b) const {
      return a.pos < b.pos;
    }
  };

  struct ssa_iterator {
    ssa_entry * e;
    ssa_entry * end;
  ssa_iterator(ssa_entry * e, ssa_entry * end) 
  : e(e), end(end) {}
    void next() {
      e++; 
    }
    int has_next() {
      return e < end;
    }
  };

  struct ssa_reverse_iterator {
    ssa_entry * c;
    ssa_entry * begin;
  ssa_reverse_iterator(ssa_entry * c, ssa_entry * begin) 
  : c(c), begin(begin) {}
    void next() {
      c--; 
    }
    int has_next() {
      return c >= begin;
    }
  };

  struct sampled_suffix_array {
    /* number of entries filled */
    idx_t n;
    /* number of elements a can hold; a power of two */
    idx_t sz;
    idx_t log_sz;			/* log(sz) */
    ssa_entry * a;
    ssa_entry * begin;		/* first entry after sort having ->pos != -1 */
    ssa_entry * end;		/* first entry after sort having ->pos != -1 */

    void init(idx_t sz_req, mallocator& mem, idx_t init_gran=10000) {
      n = 0;
      /* get a minimum power of two > sz_req */
      log_sz = calc_log_sample_sa_sz(sz_req);
      sz = (1 << (log_sz + 1));
      /* allocate the array */
      a = mem.new_<ssa_entry>(sz, mem_reason_sample_sa);
      /* parallel (doall) */
      // for(idx_t i = 0l; i < sz; i++) 
      pfor(idx_t, i, (idx_t)0, sz, init_gran) {
	a[i].rank = a[i].pos = -1;
      } end_pfor;
      begin = 0;
      end = a + sz;
    }
    void fini(mallocator& mem) {
      if (a) {
	mem.delete_(a, sz, mem_reason_sample_sa);
      }
    }

    void sort_by_pos(mallocator& mem, 
		     idx_t sort_rec_threshold=30,
		     idx_t merge_rec_threshold=1000) {
      assert(begin == 0);
      stat.start(ts_event_ssa_sort_by_pos);
      ssa_pos_lt lt;
      /* parallel sort */
      psort(a, a + sz, lt, mem, mem_reason_sort_sample_sa,
	    sort_rec_threshold, merge_rec_threshold);
      /* find first non -1 entry.
	 TODO: use binary search */
      ssa_entry * p;
      for (p = a; p < a + sz; p++) {
	if (p->pos != -1) break;
      }
      begin = p;
      end = a + sz;
      stat.end(ts_event_ssa_sort_by_pos);
    }

    /* return an iterator pointing to the largest 
       element whose ->pos <= j */
    ssa_reverse_iterator pos_reverse_begin_from(idx_t j) {
      /* TODO : binary search */
      assert(begin);
      ssa_entry * c = 0;
      for (c = end - 1; c >= begin; c--) {
	if (c->pos <= j) break;
      }
      return ssa_reverse_iterator(c, begin);
    }

    ssa_iterator pos_begin(mallocator& mem) {
      if (begin == 0) {
	sort_by_pos(mem);
      }
      return ssa_iterator(begin, end);
    }

    idx_t hash(idx_t r) {
      return r & (sz - 1);
    }

    /* return position of the rank */
    idx_t get(idx_t rank) {
      assert(begin == 0);
      stat.start(ts_event_ssa_get);
      assert(rank >= 0);
      idx_t h = hash(rank);
      for (int t = 0; t < 2; t++) {
	idx_t begin = (t == 0 ?  h : 0);
	idx_t end   = (t == 0 ? sz : h);
	for (idx_t i = begin; i < end; i++) {
	  if (a[i].rank == rank) {
	    stat.end(ts_event_ssa_get);
	    idx_t pos = a[i].pos;
	    assert(pos != -1);
	    return pos;
	  }
	  if (a[i].rank == -1) {
	    stat.end(ts_event_ssa_get);
	    return -1;
	  }
	}
      }
      stat.end(ts_event_ssa_get);
      return -1;
    }

    idx_t count() {
      idx_t c = 0;
      for (idx_t i = 0; i < sz; i++) {
	if (a[i].rank != -1) {
	  c++;
	}
      }
      return c;
    }

    int put(idx_t rank, idx_t pos) {
      assert(begin == 0);
      stat.start(ts_event_ssa_put);
      assert(rank >= 0);
      idx_t h = hash(rank);
      for (int t = 0; t < 2; t++) {
	idx_t begin = (t == 0 ?  h : 0);
	idx_t end   = (t == 0 ? sz : h);
	for (idx_t i = begin; i < end; i++) {
	  volatile bwt::ssa_entry * ai = &a[i];
	  if (ai->rank == rank) {
	    assert(ai->pos == pos);
	    stat.end(ts_event_ssa_put);
	    return 0;
	  }
	  if (ai->rank == -1) {
	    if (cmp_and_set_idx(&ai->rank, -1, rank)) {
	      ai->pos = pos;
	      assert(n < sz);
	      n++;
	      stat.end(ts_event_ssa_put);
	      return 1;
	    }
	  }
	}
      }
      check(0);			/* should not happen */
      stat.end(ts_event_ssa_put);
      return 123456;
    }
  };

}
