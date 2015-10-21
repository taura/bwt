/* 
 * conc.h --- concurrency primitives
 */
#pragma once

#include <stdint.h>
#include "typedefs.h"
#include "util.h"

#define parallel_model_none       0
#define parallel_model_native_tbb 1
#define parallel_model_task       2

#if !defined(parallel_model)
#if TO_SERIAL
#define parallel_model            parallel_model_none
#elif TO_TBB || TO_MTHREAD_NATIVE
#define parallel_model            parallel_model_task
#else
#warning "parallel_model not defined; set to parallel_model_none"
#warning "consider giving -Dparallel_model=parallel_model_{none,native_tbb,task} in the command line"
#endif
#endif

#if parallel_model == parallel_model_none
#include <algorithm>
#define decl_task_group int
#define tg_run(tg, X) X
#define tg_wait(tg)   (void)tg

#elif parallel_model == parallel_model_native_tbb
#include <tbb/task_group.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#define decl_task_group tbb::task_group
#define tg_run(tg, X) tg.run([&] { X; })
#define tg_wait(tg)   tg.wait()

#elif parallel_model == parallel_model_task
#include <common.h>
#include <mtbb/task_group.h>
#define decl_task_group mtbb::task_group
#define tg_run(tg, X) tg.run_([&] { X; }, __FILE__, __LINE__)
#define tg_wait(tg) tg.wait_(__FILE__, __LINE__)

#else
#error "define parallel_model"
#endif

#include "task_parallel_sort.h"

namespace bwt {

  bool cmp_and_set_8(volatile uint8_t * p, uint8_t x, uint8_t y) {
    return __sync_bool_compare_and_swap(p, x, y);
  }

  bool cmp_and_set_idx(volatile idx_t * p, idx_t x, idx_t y) {
    return __sync_bool_compare_and_swap(p, x, y);
  }

  void atomic_inc(volatile idx_t * p) {
    __sync_fetch_and_add(p, 1);
  }

  /* pfor macro to absorb difference between
     serial for loop, tbb::parallel_for, and
     tasks built on top of tasks;
     we may potentially include cilk in future. 

     pfor(int, i, a, b, g) {
     S
     } end_pfor;

  */

#if parallel_model == parallel_model_none
  /* 
   * pfor and pfor step for serial execution
   */

#define pfor(type, var, init, fin, gran)	\
  (void)gran;					\
  for(type var = init; var < fin; var++)
#define end_pfor

#define pfor_step(type, var, init, fin, step)	\
  for(type var = init; var < fin; var += step)
#define end_pfor_step


#elif parallel_model == parallel_model_native_tbb
  /* 
   * pfor and pfor step for native tbb
   */

#define pfor(type, var, init, fin, gran)				\
  tbb::parallel_for(tbb::blocked_range<type>(init,fin,gran), [&] (const tbb::blocked_range<type>& __range__) { for(type var = __range__.begin(); var < __range__.end(); var++) 
#define end_pfor })

#define pfor_step(type, var, init, fin, step)		\
  tbb::parallel_for(init, fin, step, [&] (type var) 
#define end_pfor_step )

#elif parallel_model == parallel_model_task
  /* 
   * pfor and pfor step with task parallel systems
   */

#define pfor(type, var, init, fin, gran)		\
  pfor_task_gran(init, fin, gran, __FILE__, __LINE__, [&] (type var) 

#define end_pfor )

#define pfor_step(type, var, init, fin, step)		\
  pfor_task_step(init, fin, step, __FILE__, __LINE__, [&] (type var) 

#define end_pfor_step )

  template<typename Index, typename Func>
    struct pfor_task_callable {
      Index first;
      Index a;
      Index b;
      Index step;
      Index gran;
      const char * file;
      int line;
      const Func & f;
    pfor_task_callable(Index first_, Index a_, Index b_, 
		       Index step_, Index gran_, 
		       const char * file, int line, 
		       const Func & f_) :
      first(first_), a(a_), b(b_), step(step_), gran(gran_), 
	file(file), line(line), f(f_) {}
      void operator() () {
	pfor_task(first, a, b, step, gran, file, line, f);
      }
    };

  template<typename Index, typename Func>
    void pfor_task(Index first, 
		   Index a, Index b, Index step, Index gran,
		   const char * file, int line, 
		   const Func& f) {
    if (b - a <= gran) {
      for (Index i = a; i < b; i++) {
	f(first + i * step);
      }
    } else {
      decl_task_group tg;
      const Index c = a + (b - a) / 2;
      tg.run_(pfor_task_callable<Index,Func>(first, a, c, step, gran, 
					     file, line, f),
	      file, line);
      pfor_task(first, c, b, step, gran, file, line, f);
      tg.wait_(file, line);
    }
  }
  
  template<typename Index, typename Func>
    void pfor_task_gran(Index first, Index last, Index gran, 
			const char * file, int line, 
			const Func& f) {
    assert(gran >= 1);
    pfor_task((Index)0, first, last, (Index)1, gran, file, line, f);
  }

  template<typename Index, typename Func>
    void pfor_task_step(Index first, Index last, Index step, 
			const char * file, int line, 
			const Func& f) {
    pfor_task(first, (Index)0, (last - first + step - 1) / step, step, (Index)1,
	      file, line, f);
  }

#else

#error "define parallel_model"

#endif

  /* ----------------------
     parallel memcpy 
     ---------------------- */
  void * pmemcpy(void * dest_, const void * src_, size_t n_, size_t gran) {
    char * dest = (char *)dest_;
    char * src = (char *)src_;
    /* todo: optimize */
    pfor(size_t, i, (size_t)0, n_, gran) {
      dest[i] = src[i];
    } end_pfor;
    return dest_;
  }

  /* prefix sum implementation */

  idx_t serial_prefix_sum(idx_t * x, idx_t a, idx_t b, idx_t t) {
    for (idx_t i = a; i < b; i++) {
      idx_t xi = x[i];
      x[i] = t;
      t += xi;
    }
    return t;
  }

#if parallel_model != parallel_model_none

  /* aux procedures for prefix sum */
  void make_sum(idx_t * x, idx_t a, idx_t b, 
		idx_t * s, idx_t p, idx_t q, idx_t g) {
    if (b - a <= g) {
      assert(q - p > 0);
      idx_t t = 0;
      for (idx_t i = a; i < b; i++) {
	t += x[i];
      }
      s[p] = t;
    } else {
      decl_task_group tg;
      idx_t n = (b - a) / 2;
      idx_t c = a + n;
      idx_t r = p + (4 * n) / g;
      tg_run(tg, make_sum(x, a, c, s, p + 1, r, g));
      make_sum(x, c, b, s, r,     q, g);
      tg_wait(tg);
      s[p] = s[p+1] + s[r];
    }
  }

  idx_t apply_sum(idx_t * x, idx_t a, idx_t b, 
		  idx_t * s, idx_t p, idx_t q, idx_t g,
		  idx_t t) {
    if (b - a <= g) {
      return serial_prefix_sum(x, a, b, t);
    } else {
      decl_task_group tg;
      idx_t n = (b - a) / 2;
      idx_t c = a + n;
      idx_t r = p + (4 * n) / g;
      idx_t t1 = t + s[p+1];
      tg_run(tg, apply_sum(x, a, c, s, p + 1, r, g, t));
      idx_t res =  apply_sum(x, c, b, s, r,     q, g, t1);
      tg_wait(tg);
      return res;
    }
  }

#endif

  idx_t prefix_sum(idx_t * x, idx_t a, idx_t b, idx_t gran,
		   mallocator& mem, mem_reason_kind_t reason) {
#if parallel_model == parallel_model_none
    (void)gran;
    (void)mem;
    (void)reason;
    return serial_prefix_sum(x, a, b, 0);
#else
    /* TODO : allocate s on stack when small? 
       there are cases make_sum is useless */
    idx_t m = max((4 * (b - a)) / gran - 1, 1);
    /* receive s rather than allocate */
    idx_t * s = mem.new_<idx_t>(m, reason);
    make_sum(x, a, b, s, 0, m, gran);
    idx_t r =  apply_sum(x, a, b, s, 0, m, gran, 0);
    mem.delete_(s, m, reason);
    return r;
#endif
  }

  /* parallel sorting */


  template<typename T, typename Compare>
    void psort(T * a_beg, T * a_end, const Compare& lt,
	       mallocator& mem, mem_reason_kind_t reason,
	       idx_t sort_rec_threshold, idx_t merge_rec_threshold) {
#if parallel_model == parallel_model_none
    (void)mem;
    (void)reason;
    (void)sort_rec_threshold;
    (void)merge_rec_threshold;
    std::sort(a_beg, a_end, lt);
#elif parallel_model == parallel_model_native_tbb
    (void)sort_rec_threshold;
    (void)merge_rec_threshold;
    tbb::parallel_sort(a_beg, a_end, lt);
#elif parallel_model == parallel_model_task
    task_parallel_sort(a_beg, a_end, lt, mem, reason,
		       sort_rec_threshold, merge_rec_threshold);
#else

#error "define parallel_model"

#endif
  }

}
