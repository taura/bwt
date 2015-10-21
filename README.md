# bwt
parallel burrows wheeler transform

## Introduction
This repository contains an implementation of
parallel and memory efficient BWT by Hayashi
and Taura:

Parallel and Memory-efficient Burrows-Wheeler
Transform.  BigData In Bioinformatics and
Health Care Informatics 2013.

The original implementation prior to
publication was not made public in a public
repository.  

Implementation in this directory is a
reworked version after the publication, and
the work contained in the directory is still
in progress.

## How to use

### For impatients

```
cd src
cp Makefile.in Makefile
make
```
(I wish I make it autotool-compliant, but currently you should manually copy Makefile.in to Makefile).

It should build two programs, simplest_example and pmbwt.


### Basics
Currently, it is provided as a set of include
files, under include/ directory.
The toplevel file to include is bwt.h.

The simplest usage is this and it is in src/simplest_example.cc

```
#include "bwt.h"

  ...
  /* input string (must end with null character) */
  const bwt::alpha_t * T = (bwt::alpha_t *)"mississippi";
  /* length must include the terminating null character */
  int n = strlen((const char *)T) + 1;
  /* allocate space for BWT */
  bwt::alpha_t * L = new bwt::alpha_t[n];

  /* set various parameters (tuning knobs) */
  bwt::bwt_opt opt;
  opt.set_defaults(T, n);
  bwt::mallocator mem(n, opt);

  /* do the real work */
  bwt::pmbwt(T, n, L, 0, mem, opt);
  ...
```

Remarks:
* Everything is defined in namespace bwt, so you should prefix all names with bwt::
* A character type is named alpha_t and it is unsigned char, so that 0 is the minimum.  So you need to explicitly cast the usual array of (signed) chars to bwt::alpha_t* (a bit awkward).
* A null character indicates the end of input characters.
* You need to pass the length of the input, include the null character (so strlen(T) + 1).
* Although the above example uses "int" as the length of input, it internally uses a 64 bit integer and it is named idx_t. 
* Knobs you can turn are defined in bwt_opt struct.

### Requirements for compilation

* Serial compilation needs a C++ compiler
* Parallel compilation using Intel TBB needs a C++ compiler that supports C++ lambda expression ([&] { ... }) and TBB library.  GNU C++ version >= 4.5 supports C++ lambda expression. 

Serial compilation:
```
cd src
g++ -Dparallel_model=parallel_model_none -I../include simplest.cc
```

Parallel compilation with Intel TBB:
```
cd src
g++ -std=c++0x -Dparallel_model=parallel_model_native_tbb -I../include simplest.cc -ltbb -ltbbmalloc
```
(depending on where you install TBB, you may need to add suitable -L and -Wl,-R options)

### Other programs you might want to compile

* src/pmbwc.cc is a command line program that generates a large random text and then build its BWT; it supports lots of parameters to play with.  It can be built by the Makefile.in in the directory (you must copy it into Makefile by yourself; sorry for the awkwardness):
```
cd src
cp Makefile.in Makefile
'''edit Makefile accordingly if you want to change options'''
make
```
You may want to run ./pmbwt --help to see what you can experiment with.
In particular, you can specify the number of cores by --n_workers.

How to specify the number of cores depends on the parallel framework you use.  Intel TBB, for example, uses all cores by default.  Pmbwt, however, changes this behavior to using whichever number of cores specified with --n_workers option, and its default is 1.  This behavior is written in pmbwt.cc and not built into the bwt library itself.  So, if you use the library from within your program, specifying the number of cores you use is up to you (or whichever underlying parallel framework you use).

* src/test/*.cc are standalone test programs for various components constituing the library.

## Discrepancy with the paper

* Currently, it does not make a suffix array of sampled characters from the whole input (as done in the original implementation) to guarantee a single string comparison does not take too long; the current implementation will thus be slow with repetitive text.  Note that it does compute sampled suffix arrays for partial strings to enable mergint two BWTs.
* It currently relies on new and delete; depending on how they internally manage memory, the memory usage may not be what you might expect (especially in parallel execution).

## Other TODOs

* Make it a template library so that you can easily change the type definition of characters and integers used to index (32 bit or 64 bit)
* Currently everything is defined in a header file, so including it from separately compiled C++ sources will result in linker errors (multiply defined symbols).  Making it template will solve the problem.
* Get rid of new/delete completely to operate within a pre-allocated memory.
* add ./configure script to streamline installation (currently you manually copy Makefile.in into Makefile)
* More performance and functionality testing

## Other Remarks

* During refactoring and reworking, I realized that doing recursions in parallel is a bit pointless, because if we do this with P processors, it increases the memory usage for temporary suffix arrays (created at leaves of the recursion) P times.  If we want to keep the total memory usage for temporary suffix arrays independent from the number of processors, we must accordingly make the leaf computation small (i.e., operate on smaller ranges).  Doing it will increase the depth of recursions and the number of times we must merge BWTs along the way, increasing the total computation so much that it may offset speedup due to parallelism.  Rather than making leaves smaller, I am leaning toward a strategy that makes leaves constant (independent from the number of processors) and only parallelize within a leaf and within a single merge operation.  I haven't fully analyzed the pros and cons of each strategy and experimental evaluations are certainly needed.
* The code is very fresh and perhaps a bit premature, so there may be things I am able to improve relatively quickly upon your request.  Please feel free to contact if you discover anything.
