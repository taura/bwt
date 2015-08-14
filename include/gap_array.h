/* 
 * gap_array.h
 */

#pragma once

#include <stdint.h>

#include "typedefs.h"
#include "util.h"
#include "conc.h"

namespace bwt {

  struct overflow_entry {
    idx_t i;			/* index */
    idx_t c;			/* count */
  };

  struct gap_array {
    idx_t a;
    idx_t b;
    idx_t m;		  /* size of overflowed array */
    uint8_t     * small;	/* counts up to 0 ... 254 */
    overflow_entry * overflow;	/* 255 or larger */
    idx_t * ps;		/* prefix sum */
    idx_t ps_a;
    idx_t ps_b;
    idx_t sum_interval;

    /* gap_array for [a:b+1] */
    void init(idx_t a_, idx_t b_, 
	      idx_t sum_interval_=128, 
	      idx_t init_gran=10000) {
      stat.start(ts_event_gap_array_init);
      a = a_;
      b = b_;
      sum_interval = sum_interval_;
      m = calc_overflow_array_sz(b - a + 1);
      ps_a = a / sum_interval;
      ps_b = b / sum_interval;

      small = new_<uint8_t>(b - a + 1, "gap small") - a;
      overflow = new_<overflow_entry>(m, "gap overflow");
      ps = new_<idx_t>(ps_b - ps_a + 1, "gap prefix sum") - ps_a;

      /* parallel (doall) */
      // for(idx_t i = a; i < b + 1; i++)
      pfor(idx_t, i, a, b + 1, init_gran) {
	small[i] = 0;
      } end_pfor;
      /* parallel (doall) */
      // for(idx_t i = 0; i < m; i++)
      pfor(idx_t, i, (idx_t)0, m, init_gran) {
	overflow[i].i = -1;
	overflow[i].c = -1;
      } end_pfor;
#if 0
      ps = 0;
      ps_a = ps_b = 0;
#endif
      stat.end(ts_event_gap_array_init);
    }

    void fini() {
      stat.start(ts_event_gap_array_fini);
      delete_(small + a, b - a + 1, "gap small");
      delete_(overflow, m, "gap overflow");
      delete_(ps + ps_a, ps_b - ps_a + 1, "gap prefix sum");
      stat.end(ts_event_gap_array_fini);
    }

    void show() {
      printf("range: [%ld,%ld] %ld elements\n", a, b, b - a + 1);
      printf("overflow: %ld elements\n", m);
      printf("prefix sum: [%ld,%ld] %ld elements\n", 
	     ps_a, ps_b, ps_b - ps_a + 1);
      printf("sum interval: %ld\n", sum_interval);
      if (small) {
	printf("small:\n");
	for (idx_t i = a; i <= b; i++) {
	  if (small[i]) printf("%6ld : %d\n", i, small[i]);
	}
      }
      if (overflow) {
	printf("overflow:\n");
	for (idx_t i = 0; i < m; i++) {
	  if (overflow[i].i != -1) {
	    printf("%6ld : %ld %ld\n", i, overflow[i].i, overflow[i].c);
	  }
	}
      }
    }

    /* smallest power of two, not smaller than n/128 */
    idx_t calc_overflow_array_sz(idx_t n) {
      idx_t y = 128;	/* 128 x */
      while (y < n) {
	y *= 2;
      }
      return y / 128;
    }

    /* try to insert to the overflow array the 
       association i -> v;
       somebody may be competing.
    */
    idx_t insert_overflow_(idx_t i) {
      idx_t h = i & (m - 1);
      /* scan the array from h to the end
	 and then beginning to h - 1 */
      for (int k = 0; k < 2; k++) {
	idx_t begin = (k == 0 ? h : 0);
	idx_t end   = (k == 0 ? m : h);
	for (idx_t j = begin; j < end; j++) {
	  volatile idx_t * p = &overflow[j].i;
	  idx_t i0 = *p;
	  if (i0 == -1) {
	    if (cmp_and_set_idx(p, i0, i)) {
	      return j;		/* I won */
	    }
	  }
	  /* either this entry is already
	     occupied or I lost in
	     cmp_and_swap; either case,
	     it should not happen the winner
	     just inserted i */
	  assert(*p != i);
	}
      }
      /* should not reach here */
      fprintf(stderr, "BUG: gap_array could not insert an entry for %ld\n", i);
      show();
      check(0);
      return 0;
    }

    /* try to insert to the overflow array the association i -> v */
    int insert_overflow(idx_t i) {
      stat.start(ts_event_gap_array_insert_overflow);
      int r = insert_overflow_(i);
      stat.end(ts_event_gap_array_insert_overflow);
      return r;
    }

    /* find an entry for i. when we find an 
       empty entry on the way, if wait_on_empty 
       is one, it waits until it becomes 
       non empty (assuming someone else may 
       be inserting i into the array);
       if wait_on_empty is zero, it immediately
       signals an error. this is used when
       you know threads who have inserted it
       have done the insertion. */
    idx_t try_find_overflow_(idx_t i, int wait_on_empty) {
      idx_t h = i & (m - 1);
      /* scan the array from h to the end
	 and then beginning to h - 1 */
      for (int k = 0; k < 2; k++) {
	idx_t begin = (k == 0 ? h : 0);
	idx_t end   = (k == 0 ? m : h);
	for (idx_t j = begin; j < end; j++) {
	  idx_t i0 = overflow[j].i;
	  if (wait_on_empty) {
	    volatile idx_t * p = &overflow[j].i;
	    for (idx_t n_tries = 0; i0 == -1 && n_tries < 100000; n_tries++) {
	      i0 = *p;
	    }
	  }
	  if (i0 == i) {
	    return j;		/* found */
	  } else if (i0 == -1) {
	    /* we find -1 before i;
	       it may be the case someone is just trying
	       to install it, but we have waited 10000 
	       times above */
	    fprintf(stderr, 
		    "BUG: gap_array could not find the entry for %ld (not full)\n", 
		    i);
	    show();
	    check(0);
	  }
	}
      }
      /* should not reach here */
      fprintf(stderr, 
	      "BUG: gap_array could not find the entry for %ld (full)\n", 
	      i);
      show();
      check(0);
      return -1;
    }

    idx_t find_overflow(idx_t i) {
      stat.start(ts_event_gap_array_find_overflow);
      idx_t j = try_find_overflow_(i, 0);
      stat.end(ts_event_gap_array_find_overflow);
      return j;
    }

    idx_t wait_overflow(idx_t i) {
      stat.start(ts_event_gap_array_wait_overflow);
      idx_t j = try_find_overflow_(i, 1);
      stat.end(ts_event_gap_array_wait_overflow);
      return j;
    }

    /* atomically increment gap[i] */
    void inc_(idx_t i) {
      assert(a <= i);
      assert(     i <= b);
      volatile uint8_t * p = &small[i];
      /* we shall not lose 256 times, as each time I lose,
	 small[i] should increase at least by one */
      for (long n_tries = 0; n_tries < 256; n_tries++) {
	uint8_t x = *p;
	if (x < 254) {
	  if (cmp_and_set_8(p, x, x + 1)) return;
	} else if (x == 254) {
	  if (cmp_and_set_8(p, x, x + 1)) {
	    /* I am the one who inserts the 
	       entry for i in the overflow array. */
	    idx_t j = insert_overflow(i);
	    assert(overflow[j].i == i);
	    assert(overflow[j].c == -1);
	    overflow[j].c = x + 1;
	    return;
	  }
	} else { 
	  assert(x == 255);
	  /* there must be somebody who set 255
	     to small[i], so, either the
	     overflow entry for i is already
	     there or somebody must be working
	     right now to install it */
	  idx_t j = wait_overflow(i);
	  assert(overflow[j].i == i);
	  volatile idx_t * q = &overflow[j].c;
	  for (idx_t m_tries = 0; m_tries < 100000; m_tries++) { 
	    if (*q != -1) {
	      atomic_inc(q);
	      return;
	    }
	  }
	  fprintf(stderr, 
		  "BUG: overflow[%ld].c remains -1 for a"
		  " long time after overflow[%ld].i == %ld (%ld)\n", 
		  j, j, i, overflow[j].i);
	  check(0);
	}
      }
      /* this is possible, but it is more likely 
	 to be a bug */
      fprintf(stderr, "BUG: incrementing small[%ld] failed 256 times\n", i);
      show();
      check(0);
    }

    void inc(idx_t i) {
      stat.start(ts_event_gap_array_inc);
      inc_(i);
      stat.end(ts_event_gap_array_inc);
    }

    idx_t get(idx_t i) {
      stat.start(ts_event_gap_array_get);
      assert(a <= i);
      assert(     i <= b);
      uint8_t x = small[i];
      if (x < 255) {
	stat.end(ts_event_gap_array_get);
	return x;
      } else {
	idx_t j = find_overflow(i);
	assert(overflow[j].i == i);
	idx_t y = overflow[j].c;
	assert(y != -1);
	stat.end(ts_event_gap_array_get);
	return y;
      }
    }

    /* set prefix sum */
    void set_prefix_sum() {
      stat.start(ts_event_gap_array_set_prefix_sum);
      /* we like to set 
	 ps[k] = get(a) + ... + get(g*k - 1)
	 (g = sum_gran) */
#if 0
      ps_a = a / sum_interval;
      ps_b = b / sum_interval;
      ps = new_<idx_t>(ps_b - ps_a + 1, "gap prefix sum") - ps_a;
#endif
      /* parallel (doall) */
      // for(idx_t k = ps_a; k < ps_b + 1; k++)
      pfor(idx_t, k, ps_a, ps_b + 1, (idx_t)1) {
	idx_t s = 0;
	for (idx_t i = max(a, k * sum_interval); 
	     i < min(b, (k + 1) * sum_interval);
	     i++) {
	  s += get(i);
	}
	/* at this point, 
	   ps[k] = get(k*g) + ... + get((k+1)*g-1) */
	ps[k] = s;
      } end_pfor;
      /* parallel (prefix sum) */
      prefix_sum(ps, ps_a, ps_b + 1, 1000);
      stat.end(ts_event_gap_array_set_prefix_sum);
    }

    /* gap[0] + ... + gap[rank] 
       (note: gap[rank] is included) */
    idx_t sum(idx_t rank) {
      stat.start(ts_event_gap_array_sum);
      assert(a <= rank);
      assert(     rank <= b);
      idx_t j = rank / sum_interval;
      assert(ps_a <= j);
      assert(        j <= ps_b);
      idx_t g = ps[j];
      idx_t i0 = rank - rank % sum_interval;
      if (i0 < a) {
	assert(g == 0);
	i0 = a;
      }
      for (idx_t i = i0; i <= rank; i++) {
	g += get(i);
      }
#if !defined(NDEBUG) || !NDEBUG
      idx_t gg = 0;
      for (idx_t i = a; i <= rank; i++) {
	gg += get(i);
      }
      assert(g == gg);
#endif
      stat.end(ts_event_gap_array_sum);
      return g;
    }

  };

}
