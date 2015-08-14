/* 
 * typedefs.h
 */

#pragma once
#include <stddef.h>

namespace bwt {
  /* alphabet_t is the type for alphabet.
     the code does not have a provision 
     for anything other than 1 byte int.
     the code uses \0 to terminate the 
     input string. */
  typedef unsigned char alpha_t;

  /* index type (integer addressing an element
     of a string or suffix array). 
     for any practical purpose, it is either 32 bit
     int or 64 bit int. it must be signed, as it
     uses -1 for a purpose */
  typedef ptrdiff_t idx_t;
  
  inline idx_t min(idx_t x, idx_t y) {
    return (x < y ? x : y);
  }
  
  inline idx_t max(idx_t x, idx_t y) {
    return (x > y ? x : y);
  }
  
}

