/* 
 * util.h --- various small stuff that do not 
 * belong to elsewhere
 */
#pragma once

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "typedefs.h"
#include "opts.h"

namespace bwt {

  /* check. it is similar to assert, but must
     should be used to checks that 
     you never want to skip, such as checking
     results. e.g.,
     
     must(x == y);
     
     it will terminate if x == y is not true */
  
  int check_(int cond, const char * exp, 
	     const char * file, int line, const char * func) {
    if (cond) {
      return 1;
    } else {
      fprintf(stderr, "error:%s:%d: check faild (%s) in %s\n",
	      file, line, exp, func);
      abort();
      return 0;
    }
  }
  
#define check(exp) check_(exp, #exp, __FILE__, __LINE__, __func__)

  /* profiling */

  typedef unsigned long long tsc_t;

  /* returns time stamp counter (x86-64 specific for now) */
  tsc_t get_tsc() {
    tsc_t u;
    asm volatile ("rdtsc;shlq $32,%%rdx;orq %%rdx,%%rax":"=a"(u)::"%rdx");
    return u;
  }
  
  /* a structure that keeps track of cumulative time
     of a specific procedure. 
     ts_record r;
     r.start();
     .. do something ..
     r.end();
     r.print();
  */
  struct ts_record {
    long n;
    tsc_t start_time;
    tsc_t total_time;
    const char * name;
    ts_record() {
      n = 0;
      start_time = 0;
      total_time = 0;
      name = 0;
    }
    /* record current time */
    void start() {
      assert(start_time == 0);
      start_time = get_tsc();
    }
    /* get current time, calc elapsed time from 
       the last start(), and accumulate the elapsed
       time */
    void end() {
      assert(start_time);
      total_time += (get_tsc() - start_time);
      start_time = 0;
      n++;
    }
    /* reset the record */
    void reset() {
      n = 0;
      start_time = 0;
      total_time = 0;
    }
    /* print
       event_name : total_time / no_of_times = avg_time
    */
    void print() {
      if (n) {
	printf("%30s : %12llu / %9ld = %14.2f\n", 
	       name, total_time, n, total_time / (double)n);
      }
    }
  };
  
  /* all the procedures we want to profile */
  typedef enum { 
    ts_event_ssa_sort_by_pos,
    ts_event_ssa_get,
    ts_event_ssa_put,
    ts_event_gap_array_init,
    ts_event_gap_array_fini,
    ts_event_gap_array_insert_overflow,
    ts_event_gap_array_wait_overflow,
    ts_event_gap_array_find_overflow,
    ts_event_gap_array_inc,
    ts_event_gap_array_get,
    ts_event_gap_array_set_prefix_sum,
    ts_event_gap_array_sum,
    ts_event_bwt_leaf,
    ts_event_bwt_sample_sa,
    ts_event_bwt_count_alphabets,
    ts_event_bwt_build_wavelet_matrix,
    ts_event_bwt_bin_search,
    ts_event_bwt_LF_map_mod,
    ts_event_sbv_rank1,
    ts_event_wm_rank,
    ts_event_bwt_LF_map_mod2,
    ts_event_bwt_sa,
    ts_event_bwt_build_gap,
    ts_event_bwt_merge_loop,
    ts_event_bwt_memcpy,
    ts_event_bwt_resample_sa,
    ts_event_bwt_merge_bwt,
    ts_event_bwt_range_rec,
    ts_event_bwt_fini,
    ts_event_pmbwt,

    n_ts_events
  } ts_event_kind_t;

  /* the array of records */
  struct ts_stat {
    int level;
    ts_record events[n_ts_events];
    ts_stat() { 
      level = 0;
#define set_event_name(x) do { events[ts_event_ ## x].name = #x; } while(0)
      set_event_name(ssa_sort_by_pos);
      set_event_name(ssa_get);
      set_event_name(ssa_put);
      set_event_name(gap_array_init);
      set_event_name(gap_array_fini);
      set_event_name(gap_array_insert_overflow);
      set_event_name(gap_array_wait_overflow);
      set_event_name(gap_array_find_overflow);
      set_event_name(gap_array_inc);
      set_event_name(gap_array_get);
      set_event_name(gap_array_set_prefix_sum);
      set_event_name(gap_array_sum);
      set_event_name(bwt_leaf);
      set_event_name(bwt_sample_sa);
      set_event_name(bwt_count_alphabets);
      set_event_name(bwt_build_wavelet_matrix);
      set_event_name(bwt_bin_search);
      set_event_name(bwt_LF_map_mod);
      set_event_name(sbv_rank1);
      set_event_name(wm_rank);
      set_event_name(bwt_LF_map_mod2);
      set_event_name(bwt_sa);
      set_event_name(bwt_build_gap);
      set_event_name(bwt_merge_loop);
      set_event_name(bwt_memcpy);
      set_event_name(bwt_resample_sa);
      set_event_name(bwt_merge_bwt);
      set_event_name(bwt_range_rec);
      set_event_name(bwt_fini);
      set_event_name(pmbwt);
#undef set_event_name
    }

    /* record timestamp for a specific procedure k */
    void start(ts_event_kind_t k) {
      if (level>=1) {
	events[k].start();
      }
    }
    /* get current time and accumulate the
       time spent in a specific procedure k */
    void end(ts_event_kind_t k) {
      if (level>=1) {
	events[k].end();
      }
    }
    /* reset all counters */
    void reset() {
      if (level>=1) {
	for (int k = 0; k < n_ts_events; k++) {
	  events[k].reset();
	}
      }
    }
    /* print all counters */
    void print() {
      if (level>=1) {
	for (int k = 0; k < n_ts_events; k++) {
	  events[k].print();
	}
      }
    }
  };

  /* the global variable to record events */
  ts_stat stat;

  /* 
     memory profiling to keep track of 
     the memory allocated and deallocated
  */
  typedef enum {
    mem_reason_leaf_sa,
    mem_reason_sort_leaf_sa,
    mem_reason_sample_sa,
    mem_reason_sort_sample_sa,
    mem_reason_count_alphabets,
    mem_reason_succinct_bitvec_superblock,
    mem_reason_succinct_bitvec_block,
    mem_reason_succinct_bitvec_prefix_sum_temp,
    mem_reason_wavelet_matrix_count_bits,
    mem_reason_wavelet_matrix_bitvec_array,
    mem_reason_wavelet_matrix_bitmaps,
    mem_reason_gap_small,
    mem_reason_gap_overflow,
    mem_reason_gap_prefix_sum,
    mem_reason_gap_prefix_sum_temp,
    mem_reason_workspace_to_merge,
    n_mem_reasons,
  } mem_reason_kind_t;

  enum {
    superblock_bits = (1 << 15),
    block_bits      = (1 << 6),
  };

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

  idx_t calc_log_sample_sa_sz(idx_t sz_req) {
    idx_t sz = 1;
    idx_t log_sz = 0;
    while (sz < sz_req) {
      sz *= 2;
      log_sz++;
    }
    return log_sz;
  }

  struct mallocator {
    /* time at which the first allocation happened */
    tsc_t started;
    /* total bytes that have been allocated so far */
    size_t allocated;
    /* total bytes that have been deallocated so far */
    size_t deallocated;
    /* max (allocated - deallocated) */
    size_t max_allocated;

    /* number of chars we are trying to compute bwt for */
    idx_t N;			
    int n_sample_sa;

    /* level of profiling. currently only 0 or 1.
       0 : not profiling
       1 : for profiling */
    int level;
    const char * reasons[n_mem_reasons];
    int n_allocated[n_mem_reasons];
    bwt_opt opt;

    mallocator(idx_t N, bwt_opt& opt) {
      level = 0;
      for (int k = 0; k < n_mem_reasons; k++) {
	n_allocated[k] = 0;
      }
      reset(N, opt);
#define set_reason_desc(x) do { reasons[mem_reason_ ## x] = #x; } while(0)
      set_reason_desc(leaf_sa);
      set_reason_desc(sort_leaf_sa);
      set_reason_desc(sample_sa);
      set_reason_desc(sort_sample_sa);
      set_reason_desc(count_alphabets);
      set_reason_desc(succinct_bitvec_superblock);
      set_reason_desc(succinct_bitvec_block);
      set_reason_desc(succinct_bitvec_prefix_sum_temp);
      set_reason_desc(wavelet_matrix_count_bits);
      set_reason_desc(wavelet_matrix_bitvec_array);
      set_reason_desc(wavelet_matrix_bitmaps);
      set_reason_desc(gap_small);
      set_reason_desc(gap_overflow);
      set_reason_desc(gap_prefix_sum);
      set_reason_desc(gap_prefix_sum_temp);
      set_reason_desc(workspace_to_merge);
#undef set_reason_desc
    }
    
    void reset(idx_t N_, bwt_opt& opt_) {
      started = 0;
      allocated = 0;
      deallocated = 0;
      max_allocated = 0;
      N = N_;
      opt = opt_;
      n_sample_sa = calc_n_sample_sa(0, N, opt_) + 1;
    }

    idx_t calc_n_sample_sa(idx_t a, idx_t b, bwt_opt& opt) {
      if (b - a <= opt.bwt_rec_threshold) {
	return 1;
      } else {
	idx_t c = (a + b) / 2;
	return 1 + calc_n_sample_sa(c, b, opt);
      }
    }

    /* record bytes are just allocated */
    template<class T>
    T * new_(size_t n, mem_reason_kind_t reason_idx) {
      T * p = new T[n];
      switch (reason_idx) {
      case mem_reason_leaf_sa:	/* leaf size */
	assert(n <= opt.bwt_rec_threshold);
	assert(n_allocated[reason_idx] < 1);
	break;
      case mem_reason_count_alphabets: /* constant */
	assert(n == opt.alpha_max + 1);
	assert(n_allocated[reason_idx] < 1);
	break;
      case mem_reason_wavelet_matrix_count_bits: /* n */
	assert(n <= ((N + opt.wavelet_matrix_sum_interval - 1) 
		     / opt.wavelet_matrix_sum_interval) * 2);
	assert(n_allocated[reason_idx] < 1);
	break;
      case mem_reason_wavelet_matrix_bitvec_array: /* constant */
	assert(n <= int_log2(opt.alpha_max));
	assert(n_allocated[reason_idx] < 1);
	break;
      case mem_reason_wavelet_matrix_bitmaps: /* n */
	assert(n <= (N / 2));
	assert(n_allocated[reason_idx] < 1);
	break;
      case mem_reason_succinct_bitvec_prefix_sum_temp:
	{
	  idx_t ns = N / superblock_bits + 1;
	  idx_t m = max((4 * ns) / opt.prefix_sum_gran - 1, 1);
	  (void)m;
	  assert(n <= m);
	  assert(n_allocated[reason_idx] < 1);
	  break;
	}
      case mem_reason_gap_small: /* n */
	assert(n <= (N / 2));
	assert(n_allocated[reason_idx] < 1);
	break;
      case mem_reason_gap_overflow: /* n/... */
	assert(n <= (N / 2) / 64);
	assert(n_allocated[reason_idx] < 1);
	break;
      case mem_reason_gap_prefix_sum: /* n/.../... */
	assert(n <= (N / 2) / opt.gap_sum_gran + 1);
	assert(n_allocated[reason_idx] < 1);
	break;
      case mem_reason_gap_prefix_sum_temp:
	{
	  idx_t ns = (N / 2) / opt.gap_sum_gran + 1;
	  idx_t m = max((4 * ns) / opt.prefix_sum_gran - 1, 1);
	  (void)m;
	  assert(n <= m);
	  assert(n_allocated[reason_idx] < 1);
	  break;
	}
      case mem_reason_workspace_to_merge: /* n */
	assert(n <= N);
	assert(n_allocated[reason_idx] < 1);
	break;
      case mem_reason_sort_leaf_sa: /*  */
	assert(n <= opt.bwt_rec_threshold);
	assert(n_allocated[reason_idx] < 1);
	break;
      case mem_reason_sort_sample_sa: /*  */
	{
	  idx_t ns = (1 << (1 + calc_log_sample_sa_sz(min(N, opt.ssa_n_samples))));
	  (void)ns;
	  assert(n <= ns);
	  assert(n_allocated[reason_idx] < 1);
	  break;
	}
      case mem_reason_sample_sa:
	{
	  idx_t ns = (1 << (1 + calc_log_sample_sa_sz(min(N, opt.ssa_n_samples))));
	  (void)ns;
	  assert(n <= ns);
	  assert(n_allocated[reason_idx] < n_sample_sa);
	  break;
	}
      case mem_reason_succinct_bitvec_superblock:
	{
	  idx_t ns = N / superblock_bits + 1;
	  (void)ns;
	  assert(n <= ns);
	  assert(n_allocated[reason_idx] < 8);
	  break;
	}
      case mem_reason_succinct_bitvec_block:
	{
	  idx_t nb = N / block_bits + 1;
	  (void)nb;
	  assert(n <= nb);
	  assert(n_allocated[reason_idx] < 8);
	  break;
	}
      default:
	check(0);
      }
      n_allocated[reason_idx]++;

      if (level>=1) {
	size_t bytes = n * sizeof(T);
	allocated += bytes;
	tsc_t t = 0;
	if (started) {
	  t = get_tsc() - started;
	} else {
	  started = get_tsc();
	}
	if (level>=2) {
	  printf("alloc/dealloc %9ld at %12llu -> %9ld bytes ( %s )\n", 
		 bytes, t, allocated - deallocated, reasons[reason_idx]);
	}
	update();
      }

      return p;
    }
    /* record bytes are just deallocated */
    template<class T>
    void delete_(T * p, size_t n, mem_reason_kind_t reason_idx) {
      assert(n_allocated[reason_idx] > 0);
      n_allocated[reason_idx]--;

      if (level>=1) {
	size_t bytes = n * sizeof(T);
	deallocated += bytes;
	tsc_t t = 0;
	if (started) {
	  t = get_tsc() - started;
	} else {
	  started = get_tsc();
	}
	if (level>=2) {
	  printf("alloc/dealloc %9ld at %12llu -> %9ld bytes ( %s )\n", 
		 -bytes, t, allocated - deallocated, reasons[reason_idx]);
	}
	update();
      }

      delete[] p;
    }
    /* update max_allocated */
    void update() {
      if (level>=1) {
	size_t a = allocated - deallocated;
	if (a > max_allocated) {
	  max_allocated = a;
	}
      }
    }
    /* print stat */
    void print() {
      if (level>=1) {
	printf("allocated     : %ld bytes\n", allocated);
	printf("deallocated   : %ld bytes\n", deallocated);
	printf("possible leak : %ld bytes\n", allocated - deallocated);
	printf("max allocated : %ld bytes\n", max_allocated);
      }
    }
  };

}
