#include <adiak.hpp>
#include <bits/getopt_core.h>
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <mpi.h>
#include <string>
#include <typeinfo>
#include <vector>

#include <caliper/cali-manager.h>
#include <caliper/cali.h>

#define IMPLEMENTATION_SOURCE "online"

enum SORT_LEVEL { SORTED, RANDOM, REVERSE, PERTURBED };
enum SCALABILITY { DEFAULT, STRONG, WEAK };
enum DATA_TYPE { INT, DOUBLE };

// Config comes with defaults
struct Config {
  long long size = 1 << 16;       // Default: 2^N or 2^16 elements
  int N = 16;                     // Default: 16
  SORT_LEVEL input_type = RANDOM; // This is the "input_type" - Default : RANDOM
  SCALABILITY scalability = DEFAULT;
  DATA_TYPE data_type = INT;
};

static void print_help(const char *program_name) {
  std::cerr
      << "Usage: " << program_name
      << " [options]\n"
         "  -n, --input-size N       Global input size (plain integer only) "
         "for 2^N. "
         "Default: 16 -> 2^16 \n"
         "  -t, --input-type TYPE_NUM    {sorted=0, random=1, reverse=2, "
         "1%perturbed=3}. "
         "Default: random\n"
         "  -s, --scaling MODE       {default=0, strong=1, weak=2}. Default: "
         "default. Just "
         "for the name so adiak can record it\n"
         "-d --data-type TYPE        {int=0, double=1}. Default: int"
         "  -h, --help               Show this help\n";
}

// input flags
// -n <size>
// -t <input_type_number>
// -s <scalability>
// -d <data_type>
// -h
int main(int argc, char *argv[]) {
  CALI_CXX_MARK_FUNCTION;

  MPI_Init(&argc, &argv);
  int rank = 0, procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  Config cfg;
  static option long_opts[] = {{"input-size", required_argument, nullptr, 'n'},
                               {"input-type", required_argument, nullptr, 't'},
                               {"scaling", required_argument, nullptr, 's'},
                               {"data-type", required_argument, nullptr, 'd'},
                               {"help", no_argument, nullptr, 'h'},
                               {nullptr, 0, nullptr, 0}};

  while (true) {
    int optindex = 0;
    int c = getopt_long(argc, argv, "n:t:s:d:h", long_opts, &optindex);

    // If there are no more input flags, then exit this loop
    if (c == -1)
      break;

    switch (c) {
    case 'n':
      cfg.size = atoi(optarg); // This the exponent for power of 2
      if (cfg.size <= 0) {
        if (rank == 0) {
          std::cerr << "Invalid input exponent, input exponent N must be "
                       "greater than 0"
                    << std::endl;
          print_help(argv[0]);
          MPI_Finalize();
        }
      } else
        cfg.size = 1 << cfg.size; // 2^cfg.size
      break;
    case 't':
      cfg.input_type = static_cast<SORT_LEVEL>(atoi(optarg));
      break;
    case 's':
      cfg.scalability = static_cast<SCALABILITY>(atoi(optarg));
      break;
    case 'd':
      cfg.data_type = static_cast<DATA_TYPE>(atoi(optarg));
      break;
    case 'h':
      print_help(argv[0]);
      MPI_Finalize();
      return 0;
    default: // Invalid argument or invalid first character
      print_help(argv[0]);
      MPI_Finalize();
      return 1;
    }
  }

  const long long N_global = cfg.size;
  const long long n_local_target =
      N_global / procs + ((rank < (N_global % procs)) ? 1 : 0);

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
                 "sample-sort"); // The name of the algorithm you are
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
  }

  // ---- From here: allocate/generate input & perform Sample Sort ----
  // Use N_global, n_local_target, cfg.input_type, cfg.trials, cfg.validate,
  // cfg.gather_output, cfg.out_prefix e.g., allocate global input on rank 0 (or
  // generate per rank), scatter, etc.
  int local_size = NULL;
  int *int_local_data = nullptr;
  double *double_local_data = nullptr;

  // ... your sorting code starts below ...
  switch (cfg.data_type) {
  case INT:
    int_local_data = new int[local_size];
    break;
  case DOUBLE:
    double_local_data = new double[local_size];
    break;
  default:
    std::cerr << "Error: Undefined data type in Config!" << std::endl;
    MPI_Finalize();
    return 1;
  }

  return 0;
}
