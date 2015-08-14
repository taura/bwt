/* 
 * pmbwt.cc
 */

#include <getopt.h>
#include "bwt.h"

/* options specific to this command */
struct pmbwt_opt {
  pmbwt_opt() {
    n_chars = 100000;
    seed = 918729723;
    alpha_min = 'a';
    alpha_max = 'z';
    stat_level = 0;
    mstat_level = 0;
    n_workers = 1;
    repeat = 1;
    help = 0;
  }
  bwt::idx_t n_chars;
  unsigned long seed;
  bwt::alpha_t alpha_min;
  bwt::alpha_t alpha_max;
  int stat_level;
  int mstat_level;
  int n_workers;
  int repeat;
  int help;
};

bwt::alpha_t * gen_input(bwt::idx_t n, 
			 bwt::alpha_t alpha0, bwt::alpha_t alpha1, 
			 unsigned long seed) {
  unsigned short rg[3] = { 
    (unsigned short)((seed >> 32) & 65535), 
    (unsigned short)((seed >> 16) & 65535), 
    (unsigned short)((seed >>  0) & 65535), 
  };
  bwt::alpha_t * T = new bwt::alpha_t[n];
  for (bwt::idx_t i = 0; i + 1 < n; i++) {
    T[i] = alpha0 + (nrand48(rg) % (alpha1 - alpha0 + 1));
  }
  T[n - 1] = 0;
  return T;
}

void random_init(bwt::alpha_t * L, bwt::idx_t n, unsigned short rg[3]) {
  for (bwt::idx_t i = 0; i < n; i++) {
    L[i] = (bwt::alpha_t)nrand48(rg);
  }
}

int check_equal(bwt::alpha_t * T, bwt::alpha_t * I, 
		bwt::idx_t a, bwt::idx_t b) {
  for (bwt::idx_t i = a; i < b; i++) {
    if (!bwt::check(T[i] == I[i])) return 0;
  }
  return 1;			/* OK */
}

int check_result(bwt::bwt& t, bwt::alpha_t * T, bwt::idx_t n) {
  /* space to get the original string back from bwt */
  bwt::alpha_t * I = new bwt::alpha_t[n];
  unsigned short rg[3] = { 918, 729, 723 };
  random_init(I, n, rg);
  printf("checking result ... "); fflush(stdout);
  t.ibwt(I);
  int r = check_equal(T, I, 0, n);
  if (r) {
    printf("OK\n");
  } else {
    printf("NG\n");
  }
  delete[] I;
  return r;
}

bwt::bwt pmbwt(bwt::alpha_t * T, bwt::alpha_t * L, bwt::idx_t n, bwt::bwt_opt& opt) {
  bwt::alpha_t * W = bwt::new_<bwt::alpha_t>(n, "workspace to merge");
  bwt::bwt t = bwt_rec(T, n, 0, n, L, W, opt);
  bwt::delete_(W, n, "workspace to merge");
  return t;
}


bwt::bwt stat_pmbwt(bwt::alpha_t * T, bwt::alpha_t * L, bwt::idx_t n, 
		    bwt::bwt_opt& opt, pmbwt_opt& opt2) {

  printf("building bwt of %ld characters\n", n);

  printf(" n_workers = %d\n", opt2.n_workers);
  printf(" repeat = %d\n", opt2.repeat);
  printf(" ssa_n_samples = %ld\n", opt.ssa_n_samples);
  printf(" bwt_rec_threshold = %ld\n", opt.bwt_rec_threshold);
  printf(" build_gap_segment_sz = %ld\n", opt.build_gap_segment_sz);
  printf(" gap_sum_gran = %ld\n", opt.gap_sum_gran);
  printf(" wavelet_matrix_sum_interval = %ld\n", 
	  opt.wavelet_matrix_sum_interval);
  printf(" n_approx_workers = %d\n", opt.n_approx_workers);

  printf(" input_alpha_min = %d\n", opt2.alpha_min);
  printf(" input_alpha_max = %d\n", opt2.alpha_max);
  printf(" algo_alpha_min = %d\n", opt.alpha_min);
  printf(" algo_alpha_max = %d\n", opt.alpha_max);

  printf(" stat_level = %d\n", opt2.stat_level);
  printf(" mstat_level = %d\n", opt2.mstat_level);
  printf(" assert_level = %d\n", opt.assert_level);
  printf(" seed = %ld\n", opt2.seed);

  bwt::stat.reset();
  bwt::stat.level = opt2.stat_level;
  bwt::mstat.reset();
  bwt::mstat.level = opt2.mstat_level;
  bwt::stat.start(bwt::ts_event_pmbwt);
  dr_start(0);
  bwt::tsc_t t0 = bwt::get_tsc();
  bwt::bwt t = pmbwt(T, L, n, opt);
  bwt::tsc_t t1 = bwt::get_tsc();
  dr_stop();
  bwt::stat.end(bwt::ts_event_pmbwt);
  printf("%llu clocks to build bwt for %ld chars\n", t1 - t0, n);
  bwt::stat.print();
  bwt::mstat.print();
  return t;
}

enum {
  opt_n_chars,
  opt_alpha_min,
  opt_alpha_max,
  opt_stat_level,
  opt_mstat_level,
  opt_repeat,
  opt_ssa_n_samples,
  opt_sort_rec_threshold,
  opt_merge_rec_threshold,
  opt_bwt_rec_threshold,
  opt_build_gap_segment_sz,
  opt_gap_sum_gran,
  opt_wavelet_matrix_sum_interval,
  opt_count_alpha_gran,
  opt_merge_gran,
  opt_memcpy_gran,
  opt_gap_init_gran,
  opt_ssa_init_gran,
  opt_wavelet_matrix_add_zero_count_gran,
  opt_n_approx_workers,
  opt_seed,
  opt_assert_level,
  n_opts,
};

struct option long_options[] = {
  {"n_workers",                 required_argument, 0, 'w' },
  {"alpha_min",                 required_argument, 0, opt_alpha_min },
  {"alpha_max",                 required_argument, 0, opt_alpha_max },
  {"stat_level",                required_argument, 0, opt_stat_level },
  {"mstat_level",               required_argument, 0, opt_mstat_level },
  {"repeat",                    required_argument, 0, 'r' },
  {"ssa_n_samples",             required_argument, 0, opt_ssa_n_samples },
  {"sort_rec_threshold",        required_argument, 0, opt_sort_rec_threshold },
  {"merge_rec_threshold",       required_argument, 0, opt_merge_rec_threshold },
  {"bwt_rec_threshold",         required_argument, 0, 't' },
  {"build_gap_segment_sz",      required_argument, 0, opt_build_gap_segment_sz },
  {"gap_sum_gran",              required_argument, 0, opt_gap_sum_gran },
  {"wavelet_matrix_sum_interval", 
   required_argument, 0, opt_wavelet_matrix_sum_interval },
  {"count_alpha_gran",          required_argument, 0, opt_count_alpha_gran },
  {"merge_gran",                required_argument, 0, opt_merge_gran },
  {"memcpy_gran",               required_argument, 0, opt_memcpy_gran },
  {"gap_init_gran",             required_argument, 0, opt_gap_init_gran },
  {"ssa_init_gran",             required_argument, 0, opt_ssa_init_gran },
  {"wavelet_matrix_add_zero_count_gran", 
                                required_argument, 0, opt_wavelet_matrix_add_zero_count_gran },
  {"n_approx_workers",          required_argument, 0, opt_n_approx_workers },
  {"seed",                      required_argument, 0, opt_seed },
  {"assert_level",              required_argument, 0, opt_assert_level },
  {"help",                      required_argument, 0, 'h' },
  {0,                           0,                 0, 0 }
};


void help() {
  for (int i = 0; long_options[i].name; i++) {
    fprintf(stderr, "--%s", long_options[i].name);
    if (long_options[i].val >= 'A') {
      fprintf(stderr, ",-%c", long_options[i].val);
    }
    switch (long_options[i].has_arg) {
    case required_argument:
      fprintf(stderr, " X");
      break;
    case optional_argument:
      fprintf(stderr, " [X]");
      break;
    default:
      break;
    }
    fprintf(stderr, "\n");
  }
}

int parse_args(int argc, char ** argv, bwt::bwt_opt& opt, pmbwt_opt& opt2) {
  (void)opt;
  assert(n_opts < 'A');

  while (1) {
    int option_index = 0;

    int c = getopt_long(argc, argv, "hw:t:r:", long_options, &option_index);
    if (c == -1) break;
    switch (c) {
    case opt_alpha_min:
      opt2.alpha_min = optarg[0];
      break;
    case opt_alpha_max:
      opt2.alpha_max = optarg[0];
      break;
    case opt_stat_level:
      opt2.stat_level = atoi(optarg);
      break;
    case opt_mstat_level:
      opt2.mstat_level = atoi(optarg);
      break;
    case 'r':
      opt2.repeat = atoi(optarg);
      break;
    case opt_ssa_n_samples:
      opt.ssa_n_samples = atol(optarg);
      break;
    case opt_sort_rec_threshold:
      opt.sort_rec_threshold = atol(optarg);
      break;
    case opt_merge_rec_threshold:
      opt.merge_rec_threshold = atol(optarg);
      break;
    case 't':
      opt.bwt_rec_threshold = atol(optarg);
      break;
    case opt_build_gap_segment_sz:
      opt.build_gap_segment_sz = atol(optarg);
      break;
    case opt_gap_sum_gran:
      opt.gap_sum_gran = atol(optarg);
      break;
    case opt_wavelet_matrix_sum_interval:
      opt.wavelet_matrix_sum_interval = atol(optarg);
      break;
    case opt_count_alpha_gran:
      opt.count_alpha_gran = atol(optarg);
      break;
    case opt_merge_gran:
      opt.merge_gran = atol(optarg);
      break;
    case opt_memcpy_gran:
      opt.memcpy_gran = atol(optarg);
      break;
    case opt_gap_init_gran:
      opt.gap_init_gran = atol(optarg);
      break;
    case opt_ssa_init_gran:
      opt.ssa_init_gran = atol(optarg);
      break;
    case opt_wavelet_matrix_add_zero_count_gran:
      opt.wavelet_matrix_add_zero_count_gran = atol(optarg);
      break;
    case opt_n_approx_workers:
      opt.n_approx_workers = atoi(optarg);
      break;
    case opt_seed:
      opt2.seed = atol(optarg);
      break;
    case opt_assert_level:
      opt.assert_level = atoi(optarg);
      break;
    case 'w':
      opt2.n_workers = atoi(optarg);
      break;
    case 'h':
      opt2.help = 1;
      return 1;
    default:
      return 0;
    }
  }

  if (optind < argc) {
    opt2.n_chars = atol(argv[optind]);
  }
  return 1;			// OK
}

#if parallel_model != parallel_model_task
#define dr_start(x) do {} while(0)
#define dr_stop()   do {} while(0)
#define dr_dump()   do {} while(0)
#endif

int main(int argc, char ** argv) {
#if parallel_model == parallel_model_task
  init_runtime(&argc, &argv);
#endif
  bwt::bwt_opt opt;
  pmbwt_opt opt2;
  if (!parse_args(argc, argv, opt, opt2)) {
    help();
    return 1;
  } else if (opt2.help) {
    help();
    return 0;
  }
#if parallel_model == parallel_model_native_tbb
  tbb::task_scheduler_init sched_init(opt2.n_workers);
#endif
  bwt::idx_t n = opt2.n_chars;
  bwt::alpha_t * T = gen_input(n, opt2.alpha_min, opt2.alpha_max, opt2.seed);
  bwt::alpha_t * L = new bwt::alpha_t[n];

  opt.set_defaults(T, n);

  int r = 1;
  for (int repeat = 0; repeat < opt2.repeat; repeat++) {
    printf("==== repeat %d ====\n", repeat);
    bwt::bwt t = stat_pmbwt(T, L, n, opt, opt2);
    r = check_result(t, T, n);
    if (!r) break;
  }
  if (r) dr_dump();

  delete[] T;
  delete[] L;
  if (r) return 0;		// OK
  else return 1;		// NG
}
