
/* T: the original string of n chars
   a: m indexes taken from T (0 <= a[i] < n)
   b: result array
   c: scratch array
 */


void dc3(alpha_t * T, idx_t n, idx_t * a, idx_t b, idx_t m) {
  idx_t n0 = (n + 2) / 3, n1 = (n + 1) / 3, n2 = n / 3, n12 = n1 + n2;
  idx_t * SA12 = new_<idx_t>(n12, "SA12");
  idx_t j = 0;
  for (idx_t i = 0; i < n; i++) {
    if (i % 3 != 0) {
      assert(j < n12);
      SA12[j++] = i;
    }
  }
  assert(j == n1);
  //******* Step 1: Sort sample suffixes ********
  // sort mod 12 suffixes by three chars
  three_char_sort(T, n, SA12, n12);
  
  // find lexicographic names of triples and
  // write them to correct places in R
  idx_t name = 0, c0 = -1, c1 = -1, c2 = -1;
  idx_t * R = new_<idx_t>(n12, "R");
  for (idx_t i = 0; i < n12; i++) {
    idx_t k = SA12[i];
    idx_t x0 = (k < n     ? T[k]   : 0);
    idx_t x1 = (k + 1 < n ? T[k+1] : 0);
    idx_t x2 = (k + 2 < n ? T[k+2] : 0);
    if (x0 != c0 || x1 != c1 || x2 != c2) {
      name++; c0 = x0; c1 = x1; c2 = x2; 
    }
    if (SA12[i] % 3 == 1) { 
      // write to R1
      assert(SA12[i] / 3 < n1);
      R[SA12[i] / 3] = name; 
    } else { 
      // write to R2
      R[SA12[i] / 3 + n1] = name; 
    }
  }
  // recurse if names are not yet unique
  if (name < n12) {
    dc3(T, n, R, SA12, n12);
    // store unique names in R using the suffix array
    for (idx_t i = 0; i < n12; i++) 
      R[SA12[i]] = i + 1;
  } else { // generate the suffix array of R directly
    for (idx_t i = 0; i < n12; i++) 
      SA12[R[i] - 1] = i;
  }
  //******* Step 2: Sort nonsample suffixes ********
  // stably sort the mod 0 suffixes from SA12 by their first character
  idx_t * R0 = new_<idx_t>(n0, "R0");
  idx_t j = 0;
  for (idx_t i = 0; i < n12; i++) {
    /* mod 0 positions. SA12[i] < n0 chooses elements originally
       at mod 1 positions */
    if (SA12[i] < n1) {
      assert(j < n0);
      R0[j++] = 3 * SA12[i];
    }
  }
  assert(j == n1);
  if (n1 < n0) {
    assert(n1 + 1 == n0);
    R0[j++] = 3 * n1;
  }
  one_char_stable_sort(T, n, R0, n0);
  //******* Step 3: Merge ********
  // merge sorted SA0 suffixes and sorted SA12 suffixes
  for (idx_t p = 0, t = 0, k = 0; k < n; k++) {
#define GetI() (SA12[t] < n1 ? SA12[t] * 3 + 1 : (SA12[t] - n1) * 3 + 2)
    idx_t i = GetI(); // pos of current offset 12 suffix
    idx_t j = SA0[p]; // pos of current offset 0 suffix
    if (SA12[t] < n1 ? // different compares for mod 1 and mod 2 suffixes
	/* mod 1 suffix */
	leq(T[i],R[SA12[t] + n1], T[j], R[j/3]) :
	/* mod 2 suffix */
	leq(T[i],T[i+1],R[SA12[t]-n0+1], T[j],T[j+1],R[j/3+n0])) {
      // suffix from SA12 is smaller
      SA[k] = i; t++;
      if (t == n02) // done --- only SA0 suffixes left
	for (k++; p < n0; p++, k++) SA[k] = SA0[p];
    } else { // suffix from SA0 is smaller
      SA[k] = j; p++;
      if (p == n0) // done --- only SA12 suffixes left
	for (k++; t < n02; t++, k++) SA[k] = GetI();
    }
  }
}

