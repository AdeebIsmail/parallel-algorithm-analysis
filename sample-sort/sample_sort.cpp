// main.cpp (lean front-half for Section 6 scripting)
// Focuses on input size (-n) and optional input type (-t); no
// seed/scaling/impl/bitshift parsing.

#include <adiak.hpp>
#include <cctype>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

#include <caliper/cali-manager.h>
#include <caliper/cali.h>

#define IMPLEMENTATION_TYPE "sample_sort"
#define IMPLEMENTATION_SOURCE "online"

struct Config {
  long long input_size = (1LL << 20); // default N = 1,048,576
  std::string input_type = "random";  // {sorted, random, reverse, 1%perturbed}
  std::string scaling = "default";
  bool validate = true;        // final correctness check
  bool gather_output = false;  // collect result on rank 0 (can be large)
  std::string out_prefix = ""; // optional: for logs/dumps
  bool verbose = false;
};

static void usage(const char *prog) {
  std::cerr
      << "Usage: " << prog
      << " [options]\n"
         "  -n, --input-size N       Global input size (plain integer only). "
         "Default: 1048576\n"
         "  -t, --input-type TYPE    {sorted, random, reverse, 1%perturbed}. "
         "Default: random\n"
         "  -s, --scaling MODE       {strong, weak}. Default: default. Just "
         "for the name so adiak can record it\n"
         "      --gather-output      Gather globally sorted array on rank 0 "
         "(memory heavy)\n"
         "  -o, --out-prefix PATH    Optional output prefix for logs/dumps\n"
         "  -v, --verbose            Verbose prints on rank 0\n"
         "  -h, --help               Show this help\n";
}

int main(int argc, char **argv) {
  CALI_CXX_MARK_FUNCTION;

  // ---- MPI first (so Adiak can see MPI context) ----
  MPI_Init(&argc, &argv);
  int rank = 0, procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  // ---- CLI: ONLY integer input-size (no bitshifts/2^k), optional input-type
  // ----
  Config cfg;
  static option long_opts[] = {{"input-size", required_argument, nullptr, 'n'},
                               {"input-type", required_argument, nullptr, 't'},
                               {"scaling", required_argument, nullptr, 's'},
                               {"gather-output", no_argument, nullptr, 1002},
                               {"out-prefix", required_argument, nullptr, 'o'},
                               {"verbose", no_argument, nullptr, 'v'},
                               {"help", no_argument, nullptr, 'h'},
                               {nullptr, 0, nullptr, 0}};

  while (true) {
    int optidx = 0;
    int c = getopt_long(argc, argv, "n:t:r:o:vh", long_opts, &optidx);
    if (c == -1)
      break;
    switch (c) {
    case 'n': {
      // plain integer only
      const char *s = optarg;
      for (const char *p = s; *p; ++p) {
        if (!std::isdigit(static_cast<unsigned char>(*p))) {
          if (rank == 0)
            std::cerr << "--input-size must be a positive integer\n";
          MPI_Finalize();
          return 2;
        }
      }
      cfg.input_size = std::atoll(s);
      break;
    }
    case 't':
      cfg.input_type = optarg;
      break;
    case 's':
      cfg.scaling = optarg;
      break;
    case 1002:
      cfg.gather_output = true;
      break;
    case 'o':
      cfg.out_prefix = optarg;
      break;
    case 'v':
      cfg.verbose = true;
      break;
    case 'h':
      usage(argv[0]);
      MPI_Finalize();
      return 0;
    default:
      usage(argv[0]);
      MPI_Finalize();
      return 1;
    }
  }

  // Normalize & check fields
  auto lower = [](std::string s) {
    for (auto &c : s)
      c = std::tolower((unsigned char)c);
    return s;
  };
  cfg.input_type = lower(cfg.input_type);
  cfg.scaling = lower(cfg.scaling);

  if (cfg.input_type == "reverse_sorted")
    cfg.input_type = "reverse";
  if (cfg.input_type == "1percent" || cfg.input_type == "1pct" ||
      cfg.input_type == "perturbed")
    cfg.input_type = "1%perturbed";

  if (cfg.scaling != "strong" && cfg.scaling != "weak" &&
      cfg.scaling != "default") {
    std::cerr << "Unknown --scaling: " << cfg.scaling << "\n";
    MPI_Finalize();
    return 4;
  }

  if (cfg.input_type != "sorted" && cfg.input_type != "random" &&
      cfg.input_type != "reverse" && cfg.input_type != "1%perturbed") {
    if (rank == 0)
      std::cerr << "Unknown --input-type: " << cfg.input_type << "\n";
    MPI_Finalize();
    return 3;
  }
  if (cfg.input_size <= 0) {
    if (rank == 0)
      std::cerr << "--input-size must be > 0\n";
    MPI_Finalize();
    return 4;
  }

  // Derived sizes for scatter (handles imbalance via +1 on the first remainder
  // ranks)
  const long long N_global = cfg.input_size;
  const long long n_local_target =
      N_global / procs + ((rank < (N_global % procs)) ? 1 : 0);

  // ---- Adiak metadata (Section 5) ----
  adiak::init(NULL);
  adiak::launchdate();  // launch date of the job
  adiak::libraries();   // Libraries used
  adiak::cmdline();     // Command line used to launch the job
  adiak::clustername(); // Name of the cluster
  adiak::value("group_num",
               6); // The number of your group (integer, e.g., 1, 10)
  adiak::value("implementation_source",
               IMPLEMENTATION_SOURCE); // Where you got the source code of your
                                       // algorithm. choices: ("online", "ai",
                                       // "handwritten").
  if (rank == 0) {
    adiak::value("algorithm",
                 IMPLEMENTATION_TYPE); // The name of the algorithm you are
                                       // using (e.g., "merge", "bitonic")
    adiak::value("programming_model", "mpi");
    adiak::value(
        "input_type",
        cfg.input_type); // For sorting, this would be choices: ("Sorted",
                         // "ReverseSorted", "Random", "1_perc_perturbed")
    adiak::value("input_size",
                 N_global); // The number of elements in input dataset (1000)
    adiak::value("scalability",
                 scalability); // The scalability of your algorithm. choices:
                               // ("strong", "weak")
    adiak::value("num_procs", procs); // The number of processors (MPI ranks)
    adiak::value("trials", cfg.trials);
    adiak::value("gather_output", cfg.gather_output ? 1 : 0);
  }

  // ---- Optional banner ----
  if (rank == 0 && cfg.verbose) {
    std::cout << "[cfg] input_type=" << cfg.input_type
              << " trials=" << cfg.trials << " procs=" << procs
              << " N_global=" << N_global
              << (cfg.gather_output ? " gather_output=on"
                                    : " gather_output=off")
              << (cfg.validate ? " validate=on" : " validate=off")
              << (cfg.out_prefix.empty() ? ""
                                         : (" out_prefix=" + cfg.out_prefix))
              << std::endl;
  }

  // ---- From here: allocate/generate input & perform Sample Sort ----
  // Use N_global, n_local_target, cfg.input_type, cfg.trials, cfg.validate,
  // cfg.gather_output, cfg.out_prefix e.g., allocate global input on rank 0 (or
  // generate per rank), scatter, etc.
  vector<double> double_vector;
  vector<int> int_vector;

  // ... your sorting code starts below ...
}
