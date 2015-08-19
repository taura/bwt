#pragma once

#include <stdio.h>

#include "typedefs.h"

namespace bwt {

  struct bwt_opt {
    bwt_opt() {
      mem_budget = 0;
      alpha_min = 0;
      alpha_max = 0;
      ssa_n_samples = 0;
      sort_rec_threshold = 0;
      merge_rec_threshold = 0;
      bwt_rec_threshold = 0;
      build_gap_segment_sz = 0;
      count_alpha_gran = 0;
      merge_gran = 0;
      memcpy_gran = 0;
      gap_init_gran = 0;
      ssa_init_gran = 0;
      sa_init_gran = 0;
      ssa_sample_gran = 0;
      gap_sum_gran = 0;
      sa_to_bwt_gran = 0;
      wavelet_matrix_sum_interval = 0;
      wavelet_matrix_add_zero_count_gran = 0;
      wavelet_matrix_transpose_gran = 0;
      prefix_sum_gran = 0;
      n_approx_workers = 0;
      assert_level = 0;
    }

    double mem_budget;

    /* minimum alphet. if 0, calculated from the input string */
    alpha_t alpha_min;
    /* minimum alphet. if 0, calculated from the input string */
    alpha_t alpha_max;
    /* algorithm tuning parameters */

    /* 
       number of samples taken for sampled
       suffix array.  as a rule of thumb, we
       want to make its size smaller than the
       input array size. i.e., if the input has
       n characters and the size of an index is
       8 bytes, we want to make it smaller than
       n/8. if you make it larger, merging
       becomes faster.  if 0, it is set to n/32.
    */
    idx_t ssa_n_samples;

    idx_t sort_rec_threshold;
    idx_t merge_rec_threshold;

    /* 
       the number of characters below which we
       do not divide further when building
       bwt. strings shorter than this is
       trasformed directly via the usual way
       (building suffix array and then
       converting it to bwt) you do not want to
       make it too small, as doing so results
       in deep recursions, which in turn
       increase the number of times each
       character experiences LF mapping on the
       way, which makes it entire computation
       more expensive. on the other hand you do
       not make it too large, as doing so
       requires a larger memory for the suffix
       array.  if 0, or if it is too small to
       be useful, the algorithm computes it
       based on n and n_workers.  if n_approx_workers
       is also not given, it is set to a
       reasonable constant 
    */

    idx_t bwt_rec_threshold;

    /* 
       the number of bytes we below which we do
       not parallelize building gap array,
       which is a process in mergint two bwts.
       when building a gap array, we divide the
       right bwt into many segments, so as to
       make each segment shorter than
       build_gap_gran.  the position of
       positionally the last character in a
       segment (which is the character whose
       position in the left array is determined
       FIRST, as we scan backward) is found
       by a binary search. positions of all
       other characters are found by
       LF-mapping. since the binary search is
       more costly than a LF-mapping, making
       this number too small makes binary search
       more frequent.  if we make this number
       too large, you may lose parallelism.
       theoretically, we want to avoid order n
       serial computation, so we want to make
       it a constant not depending on n, or
       something like log(n).
       if 0, the algorithm computes an appropriate
       choice based on n and n_approx_workers.
       if n_approx_workers is not given (0),
       it chooses a reasonable constant. 
    */

    idx_t build_gap_segment_sz;

    /* 
       the number of bytes between which the
       prefix sum of gap array is stored; after
       a gap array is built, we know gap[i] is
       the number of characters to insert
       between i-th and (i+1)-th elements of
       the left array. we like to transform
       the gap info into the prefix sum, so 
       as to know the position to which each
       character in the right array goes.
       we of course do not have the luxury
       of having the array having as many elements
       as the original left array, so we 
       sparsely store samples of the resulting
       prefix sum array. gap_sum_gran is
       the interval at which these samples are
       taken.  if you increase this number,
       we use less memory for the prefix array,
       but determining the position of each
       character takes more time
       (roughly proportional to gap_sum_gran).
       as a rule of thum, you want to make
       8n/gap_sum_gran much smaller than n,
       so you want to make gap_sum_gran larger
       enough than 8, usually 128 or something.
       if 0, a reasonable constant is chosen
       by the algorithm.
    */

    idx_t gap_sum_gran;

    idx_t sa_to_bwt_gran;

    /* 
       when buidling a wavelet matrix of a
       string, we repeatedly partition the
       characters by a bit of each character;
       all characters having 0 in their MSBs go
       left of the array and all having 1 go
       right.  we next partition them by the
       second MSB, and the process goes on.  to
       do this, we desire to calculate, for
       each index i, the number of characters
       in s[0:i] having a 0 or 1 at the
       specified bit position.  literally
       computing this for every position is
       simply too expensive and not necessary,
       so we do this at every
       wavelet_matrix_segment_sz bytes.
    */
    idx_t wavelet_matrix_sum_interval;

    idx_t count_alpha_gran;
    idx_t merge_gran;
    idx_t memcpy_gran;
    idx_t gap_init_gran;
    idx_t ssa_init_gran;
    idx_t sa_init_gran;
    idx_t ssa_sample_gran;
    idx_t wavelet_matrix_add_zero_count_gran;
    idx_t wavelet_matrix_transpose_gran;
    idx_t prefix_sum_gran;

    /* 
       the number of workers you use, or 0 if
       you have no idea. this value is used
       only for choosing parameters for good
       performance/memory usage; the
       correctness of the algorithm never
       depends on this value.
    */

    int n_approx_workers;

    /* 
       assertion level.

       there are very expensive assertions you
       don't want to perform unless you are
       diagnosing bugs with small inputs.  when
       building with assertion on (i.e., NDEBUG
       is either undefined or set to zero), such
       expensive assertion can be turned on/off
       dynamically without recompilation.  set it
       to a value >= 2 to turn on those expensive
       assertion checks.
    */

    int assert_level;

    idx_t max_num_ssa(idx_t n, idx_t th) {
      if (n <= th) return 0;
      else {
	return 1 + max_num_ssa(n - n / 2, th);
      }
    }

    void set_defaults(const alpha_t * T, idx_t n) {
      (void)T;			/* not used at this point */
      if (alpha_min == 0 && alpha_max == 0) {
	alpha_min = 1;
	alpha_max = 255;
      }
      if (n_approx_workers == 0) {
	n_approx_workers = 256;
      }
      if (sort_rec_threshold == 0) {
	sort_rec_threshold = 30;
      }
      if (merge_rec_threshold == 0) {
	merge_rec_threshold = 10000;
      }
      if (build_gap_segment_sz == 0) {
	build_gap_segment_sz = 1024;
      }
      if (gap_sum_gran == 0) {
	gap_sum_gran = 128;
      }
      if (sa_to_bwt_gran == 0) {
	sa_to_bwt_gran = 10000;
      }
      if (wavelet_matrix_sum_interval == 0) {
	wavelet_matrix_sum_interval = 4096;
      }
      if (count_alpha_gran == 0) {
	count_alpha_gran = 10000;
      }
      if (merge_gran == 0) {
	merge_gran = 1000;
      }
      if (memcpy_gran == 0) {
	memcpy_gran = 10000;
      }
      if (gap_init_gran == 0) {
	gap_init_gran = 10000;
      }
      if (ssa_init_gran == 0) {
	ssa_init_gran = 10000;
      }
      if (sa_init_gran == 0) {
	sa_init_gran = 10000;
      }
      if (ssa_sample_gran == 0) {
	ssa_sample_gran = 10000;
      }
      if (wavelet_matrix_add_zero_count_gran == 0) {
	wavelet_matrix_add_zero_count_gran = 10000;
      }
      if (wavelet_matrix_transpose_gran == 0) {
	wavelet_matrix_transpose_gran = 10000;
      }
      if (prefix_sum_gran == 0) {
	prefix_sum_gran = 1000;
      }
      if (mem_budget == 0) {
	mem_budget = 4.0;
      }
      if (mem_budget < 2.5) {
	printf("warning: mem_budget %f too small. defaults to 4.0\n", mem_budget);
	mem_budget = 4.0;
      }
      if (mem_budget >= 2 + sizeof(idx_t) * 2) {
	bwt_rec_threshold = n;
	ssa_n_samples = 2;	/* should not matter */
      } else {
	/* allocate 80% of the remaining budget to leaf suffix array and its sorting */
	double p = 0.8;
	bwt_rec_threshold = ((mem_budget - 2) * p * n) / (sizeof(idx_t) * 2);
	bwt_rec_threshold = max(1000, bwt_rec_threshold);
	idx_t ssa_total_bytes = (mem_budget - 2) * (1 - p) * n;
	idx_t ssa_bytes = ssa_total_bytes / (1 + max_num_ssa(n, bwt_rec_threshold));
	ssa_n_samples = ssa_bytes / (sizeof(idx_t) * 2);
	ssa_n_samples = max(1000, ssa_n_samples);
      }
    }

  };

}
