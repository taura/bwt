/* 
 * succinct bit vector
 */
#pragma once

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#include "typedefs.h"
#include "util.h"
#include "conc.h"

namespace bwt {

  /*
  
    succinct_bit_vector s(bits, n);
  
    s.rank(x) returns the number of bits 
    in the first x bits of a

  */

  enum {
    x_superblock_bits = (1 << 15),
    x_block_bits      = (1 << 9),
    superblock_bits = (1 << 15),
    block_bits      = (1 << 6),
  };

  struct succinct_bit_vector {

  private:
    uint8_t * bits;
    idx_t begin;			/* bits + begin ... bits + end - 1 */
    idx_t end;
    idx_t n_superblocks;
    idx_t    * superblocks;
    idx_t n_blocks;
    uint16_t * blocks;

  public:
    succinct_bit_vector() {
      bits = 0;
      superblocks = 0;
      blocks = 0;
    }
    idx_t n_zeros;

  private:
    /* naively count the number of bits set in an uint64.
       later replace it by pop count */
    idx_t count_ones_in_uint64_portable(uint64_t u) {
      idx_t s = 0;
      for (int i = 0; i < 64; i++) {
	s += ((u >> i) & 1);
      }
      return s;
    }

    /* count the number of bits set in an uint64,
       using popcount instruction */
    idx_t count_ones_in_uint64_popcount(uint64_t u) {
      return __builtin_popcountl(u);
    }

    /* count ones in a single 64 bit word
       TODO : ifdef based on the availability of popcount */
    idx_t count_ones_in_uint64(uint64_t u) {
      return count_ones_in_uint64_popcount(u);
      // return count_ones_in_uint64_portable(u);
    }

    /* count ones in address [b,e) for aligned addresses;
       both b and e are aligned to 8 bytes (64 bits) */
    __attribute__((optimize("unroll-loops")))
    idx_t count_ones_in_range_64(uint64_t * b, uint64_t * e) {
      idx_t s = 0;
      assert((intptr_t)b % 8 == 0);
      assert((intptr_t)e % 8 == 0);
      for (uint64_t * p = b; p < e; p++) {
	s += count_ones_in_uint64(*p);
      }
      return s;
    }

    /* count ones in address [b,e) for arbitrary 
       (possibly unaligned) addresses b and e

       b64 ...... b ............. e64 ...... e
       |<---------    x   ------->|
       |<-  y  ->|                |<-  z  ->|

       we compute x + z - y  

       b ....... b64 ............ e64 ...... e
       |<----  x   ---->|
       |<-  y  ->|                |<-  z  ->|


    */
    __attribute__((optimize("unroll-loops")))
    idx_t count_ones_in_range_8(uint8_t * b, uint8_t * e) {
      /* 64 bit (8 bytes)-aligned addresses b64 <= b and e64 <= e */
      uint8_t * b64 = (uint8_t *)((intptr_t)(b + 7) - (intptr_t)(b + 7) % 8);
      uint8_t * e64 = (uint8_t *)((intptr_t)e - (intptr_t)e % 8);
      idx_t rb = 0;
      idx_t re = 0;
      idx_t s = 0;
      assert(b <= b64);
      assert(b64 - b < 8);
      assert(e64 <= e);
      assert(e - e64 < 8);

      if (b64 <= e64) {
	/* count ones in [b,b64] */
	assert(b64 <= e);
	for (uint8_t * p = b; p < b64; p++) {
	  rb = (rb << 8) + *p;
	}
	s += count_ones_in_uint64(rb);
	/* count ones in [b64,e64] */
	s += count_ones_in_range_64((uint64_t *)b64, (uint64_t *)e64);
	/* count ones in [e64,e] */
	for (uint8_t * p = e64; p < e; p++) {
	  re = (re << 8) + *p;
	}
	s += count_ones_in_uint64(re);
      } else {
	for (uint8_t * p = b; p < e; p++) {
	  rb = (rb << 8) + *p;
	}
	s += count_ones_in_uint64(rb);
      }
      return s;
    }

    /* count the number of ones in arbitrary bit positions b and e
       (note that b and e are given in the number of BITS, not bytes).

       b  ...... b8 ............. e8  ...... e
       |<---    x   --->|
       |<-  y  ->|                |<-  z  ->|

       we compute y + x + z

       we maintain we never touch bytes not containing
       the bit range [b:e)

       note that when b and e are very close, the
       point to two bit positions in a single byte.
       in that case, it might happen that e8 < b8
       (e8 + 8 == b8, to be specific)
    */
    idx_t count_ones_in_range(idx_t b_, idx_t e_) {
      assert(b_ <= e_);
      uint8_t * h = bits;
      /* bit position from the first bit of the BITS */
      idx_t b = b_ + begin;
      idx_t e = e_ + begin;
      /* get 8 bit-aligned bit positions b <= b8 and e8 <= e  */
      idx_t b8 = (b + 7) - (b + 7) % 8;
      idx_t e8 = e - e % 8;
      /* handle the case when b and e point to the same byte */
      if (e8 < b8) {
	assert(b < b8);
	assert(e8 < e);
	assert(e8 + 8 == b8);
	/* 
	   e8  b  e
	   ........
	   00000100  1 << (e - e8)
	   00100000  1 << (b - e8)
	   -) 00111000  
	*/
	idx_t c = h[b8 / 8 - 1];
	return count_ones_in_uint64(c & ((1 << (e - e8)) - (1 << (b - e8))));
      } else {
	idx_t s = 0;
	if (b < b8) 		/* y part, if there is any overlap */
	  s += count_ones_in_uint64(h[b8 / 8 - 1] >> (b + 8 - b8));
	/* x part */
	s += count_ones_in_range_8(h + b8 / 8, h + e8 / 8);
	if (e8 < e) 		/* z part, if there is any overlap */
	  s += count_ones_in_uint64(h[e8 / 8] & ((1 << (e - e8)) - 1));
	return s;
      }
    }

    /* only for test. 1 if idx-th bit is set */
    idx_t is_bit_set(idx_t idx) {
      uint8_t * h = bits;
      idx_t byte_idx = (idx + begin) / 8;
      idx_t bit_idx  = (idx + begin) % 8;
      return (h[byte_idx] >> bit_idx) & 1;
    }

  public:
    void init(uint8_t * bits_, idx_t begin_, idx_t end_) {
      /* we don't like a block not aligned to 
	 64 bit word boundaries ... */
      assert(block_bits % 64 == 0);
      /* nor a superblock not aligned to block boundaries ... */
      assert(superblock_bits % block_bits == 0);
      /* we must have the number of ones within 
	 a single superblock representable by a 
	 single 16 bit word (what we use to count 
	 ones within a single superblock) */
      assert(superblock_bits < (1UL << (sizeof(blocks[0]) * 8)));
      bits = bits_;
      begin = begin_;
      end = end_;
      idx_t n = end - begin;
      /* total number of superblocks and blocks */
      idx_t ns = n / superblock_bits + 1;
      idx_t nb = n /      block_bits + 1;
      idx_t n_blocks_per_superblock = superblock_bits / block_bits;
      n_superblocks = ns;
      n_blocks      = nb;
      superblocks = new_<idx_t>(ns, "succinct bit superblock");
      blocks      = new_<uint16_t>(nb, "succinct bit block");
      /* the first pass. count ones in each superblock
	 and each block. then we make the block counts
	 a prefix sum within each superblock */
    
      /* parallel (doall) */
      // for (idx_t s = 0; s < ns; s++) 
      pfor(idx_t, s, (idx_t)0, ns, (idx_t)1) {
	idx_t b_start = s * n_blocks_per_superblock;
	idx_t b_end   = min(b_start + n_blocks_per_superblock, nb);
	superblocks[s] = 0;
	for (idx_t b = b_start; b < b_end; b++) {
	  idx_t i_start = b * block_bits;
	  idx_t i_end   = min(i_start + block_bits, n);
	  idx_t n1 = count_ones_in_range(i_start, i_end);
	  blocks[b] = n1;
	  superblocks[s] += n1;
	}
	/* make blocks the prefix sum within this superblock */
	idx_t t = 0;
	for (idx_t b = b_start; b < b_end; b++) {
	  idx_t c = blocks[b];
	  blocks[b] = t;
	  t += c;
	}
      } end_pfor;
      /* the second path to get the prefix sum */
      /* parallel (prefix sum) */
      prefix_sum(superblocks, 0, ns, 1000);
      /* remember number of zeros for wavelet matrix */
      n_zeros = rank0(n);
    }

    void fini() {
      if (superblocks) 
	delete_(superblocks, n_superblocks, "succinct bit superblock");
      if (blocks) 
	delete_(blocks, n_blocks, "succinct bit block");
    }

    /* the number of bits in a[0:n] */
    idx_t rank1(idx_t n) {
      stat.start(ts_event_sbv_rank1);
      assert(0 <= n);
      assert(     n <= end - begin);
      assert(n / superblock_bits < n_superblocks);
      assert(n / block_bits      < n_blocks);
      idx_t c0 = superblocks[n / superblock_bits];
      idx_t c1 = blocks[n / block_bits];
      idx_t c2 = count_ones_in_range(n - n % block_bits, n);
      idx_t c012 = c0 + c1 + c2;
      assert(0 <= c0);
      assert(     c0 <= n);
      assert(0 <= c1);
      assert(     c1 <= n);
      assert(0 <= c2);
      assert(     c2 <= n);
      assert(c012 <= n);
      stat.end(ts_event_sbv_rank1);
      return c012;
    }

    /* the number of zero bits in a[0:n] */
    idx_t rank0(idx_t n) {
      return n - rank1(n);
    }

    /* only for test. supernaive way to count set bits in bits[b:e] */
    idx_t rank_slow(idx_t b, idx_t e) {
      idx_t s = 0;
      for (idx_t i = b; i < e; i++) {
	s += is_bit_set(i);
      }
      return s;
    }

  };

}
