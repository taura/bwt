/* 
 * util.h --- various small stuff that do not 
 * belong to elsewhere
 */
#pragma once

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <new>

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
  ts_record() : n(0), start_time(0), total_time(0), name(0) { }
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
  ts_stat() : level(0) { 
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
  }
    ts_record events[n_ts_events];

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
  struct mem_stat {
    /* time at which the first allocation happened */
    tsc_t started;
    /* total bytes that have been allocated so far */
    size_t allocated;
    /* total bytes that have been deallocated so far */
    size_t deallocated;
    /* max (allocated - deallocated) */
    size_t max_allocated;
    /* level of profiling. currently only 0 or 1.
       0 : not profiling
       1 : for profiling */
    int level;

  mem_stat() : started(0), allocated(0), deallocated(0), 
      max_allocated(0), level(0) {}
    void reset() {
      started = 0;
      allocated = 0;
      deallocated = 0;
      max_allocated = 0;
    }
    /* record bytes are just allocated */
    void new_(size_t bytes, const char * reason) {
      if (level>=1) {
	allocated += bytes;
	tsc_t t = 0;
	if (started) {
	  t = get_tsc() - started;
	} else {
	  started = get_tsc();
	}
	if (level>=2) {
	  printf("alloc/dealloc %9ld at %12llu -> %9ld bytes ( %s )\n", 
		 bytes, t, allocated - deallocated, reason);
	}
	update();
      }
    }
    /* record bytes are just deallocated */
    void delete_(size_t bytes, const char * reason) {
      if (level>=1) {
	deallocated += bytes;
	tsc_t t = 0;
	if (started) {
	  t = get_tsc() - started;
	} else {
	  started = get_tsc();
	}
	if (level>=2) {
	  printf("alloc/dealloc %9ld at %12llu -> %9ld bytes ( %s )\n", 
		 -bytes, t, allocated - deallocated, reason);
	}
	update();
      }
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

  /* global variable for memory profiling */
  mem_stat mstat;

  /* wrappers for new[] and delete[] */

  /*     T * t = new T[n];
	 ==> T * t = new_<T>(n, "reason");
  */
  template<class T>
    T * new_(size_t n, const char * reason) {
    mstat.new_(sizeof(T) * n, reason);
    return new T[n];
  }

  /*     delete[] t;
	 ==> delete(t, n, "reason")
  */
  template<class T>
    void delete_(T * p, size_t n, const char * reason) {
    mstat.delete_(sizeof(T) * n, reason);
    delete[] p;
  }

}
