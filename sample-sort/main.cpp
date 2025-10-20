#include <adiak.hpp>
#include <bits/getopt_core.h>
#include <cctype>
#include <climits>
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

struct Dist {
  long long *counts;
  long long *displs;
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

// Note: the count and displ are not in bytes, but in number of elements
static void get_count_displacement(Dist &dist, long long N_global_size,
                                   int num_procs, int rank) {
  dist.counts = new long long[num_procs];
  dist.displs = new long long[num_procs];

  long long base_chunk_size = N_global_size / num_procs;
  int remainder = N_global_size % num_procs;
  long long running = 0;

  for (int p = 0; p < num_procs; p++) {
    long long c_p = base_chunk_size + (p < remainder ? 1 : 0);
    dist.counts[p] = c_p;
    dist.displs[p] = running;
    running += c_p;
  }
}

static void local_sort_mergesort(const void *local_array,
                                 const long long local_size,
                                 const DATA_TYPE d_type);

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

    // c -> Flag being parsed
    switch (c) {
    case 'n':
      cfg.N = atoi(optarg); // This the exponent for power of 2
      if (cfg.size <= 0) {
        if (rank == 0) {
          std::cerr << "Invalid input exponent, input exponent N must be "
                       "greater than 0"
                    << std::endl;
          print_help(argv[0]);
          MPI_Finalize();
        }
      } else
        cfg.size = 1 << cfg.N; // 2^N
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

  std::string scalability;
  switch (cfg.scalability) {
  case DEFAULT:
    scalability = "default";
    break;
  case WEAK:
    scalability = "weak";
    break;
  case STRONG:
    scalability = "strong";
    break;
  default:
    scalability = "default";
    break;
  }

  std::string input_type;
  switch (cfg.input_type) {
  case RANDOM:
    input_type = "Random";
    break;
  case SORTED:
    input_type = "Sorted";
    break;
  case REVERSE:
    input_type = "ReverseSorted";
    break;
  case PERTURBED:
    input_type = "1_perc_perturbed";
    break;
  }

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
    adiak::value("input_type",
                 input_type); // For sorting, this would be choices: ("Sorted",
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
  //
  // # ---------- counts & displacements (elements, not bytes) ----------
  int *int_local_data = nullptr;
  double *double_local_data = nullptr;
  Dist dist;

  // Getting everyone's local_size and their displacements in the global array
  get_count_displacement(dist, N_global, procs, rank);

  long long local_size = dist.counts[rank];
  long long global_offset = dist.displs[rank];

  //  # ---------- Distributed initialization (no root global buffer) ----------
  CALI_MARK_BEGIN("data_init_runtime");
  switch (cfg.data_type) {
  case INT:
    int_local_data = new int[local_size];

    switch (cfg.input_type) {
    case RANDOM:
      for (int i = 0; i < local_size; i++) {
        int_local_data[i] = rand() % 100000;
      }
      break;
    case SORTED:
      for (int i = 0; i < local_size; i++) {
        int_local_data[i] = static_cast<int>(global_offset + i);
      }
      break;
    case REVERSE:
      for (int i = 0; i < local_size; i++) {
        const long long g = global_offset + i;
        int_local_data[i] = static_cast<int>((N_global - 1) - g);
      }
      break;
    case PERTURBED:
      for (int i = 0; i < local_size; i++) {
        int_local_data[i] = static_cast<int>(global_offset + i);
      }

      const int num_swaps = max(1, floor(0.01 * local_size));
      for (int i = 0; i < num_swaps; i++) {
        int idx = (rand() % local_size) - 1;
        int temp = int_local_data[i];
        int_local_data[idx] = data[idx + 1];
        int_local_data[idx + 1] = temp;
      }
      break;
    }
    break;
  case DOUBLE:
    double_local_data = new double[local_size];

    switch (cfg.input_type) {
    case RANDOM:
      for (int i = 0; i < local_size; i++) {
        double_local_data[i] = rand() % 100000;
      }
      break;
    case SORTED:
      for (int i = 0; i < local_size; i++) {
        double_local_data[i] = static_cast<double>(global_offset + i);
      }
      break;
    case REVERSE:
      for (int i = 0; i < local_size; i++) {
        const long long g = global_offset + i;
        double_local_data[i] = static_cast<double>((N_global - 1) - g);
      }
      break;
    case PERTURBED:
      for (int i = 0; i < local_size; i++) {
        double_local_data[i] = static_cast<double>(global_offset + i);
      }

      const int num_swaps = max(1, floor(0.01 * local_size));
      for (int i = 0; i < num_swaps; i++) {
        int idx = (rand() % local_size) - 1;
        double temp = int_local_data[i];
        double_local_data[idx] = data[idx + 1];
        double_local_data[idx + 1] = temp;
      }
      break;
    }
    break;
  default:
    std::cerr << "Error: Undefined data type in Config!" << std::endl;
    MPI_Finalize();
    return 1;
  }
  CALI_MARK_END("data_init_runtime");

  // # ===================== Sample Sort proper =====================
  // 1) Local sort
  CALI_MARK_BEGIN("comp");
  CALI_MARK_BEGIN("comp_large");
  switch (cfg.data_type) {
  case INT: {
    local_sort_mergesort(int_local_data, local_size, cfg.data_type);

    // 2) Local sampling: pick s=P evenly spaced samples
    const size_t samples = procs;
    int *int_samples = new int[samples];
    for (int i = 0; i < samples; i++) {
      // TODO: Unsure how this can go wrong though...
      size_t idx =
          (size_t)((i + 1.0) * (local.size()) / (s + 1.0)); // evenly spaced
      int_samples[i] = int_local_data[idx];
    }

    // 3) Allgather samples (small)
    const size_t all_samples_size = samples * procs;
    int *all_samples = new int[all_samples_size];
    MPI_Allgather(int_samples, samples, MPI_INT, all_samples, samples, MPI_INT,
                  MPI_COMM_WORLD); // TODO: COM_WORLD might be wrong here...

    // 4) Choose global splitters (small)
    local_sort_mergesort(all_samples, all_samples_size, cfg.data_type);
    int *splitters = new int[procs - 1];
    for (int i = 1; i < procs; i++) {
      size_t idx = (size_t)((i * all_samples_size) / procs);
      splitters[i - 1] = all_samples[idx];
    }

    // 5) Partition LOCAL into P buckets by splitters
    int *buckets = new int[procs];
    for (int bucket_num = 0; bucket_num < procs; bucket_num++) {
      int lower = (bucket_num == 0) ? INT_MIN : splitters[bucket_num - 1];
      int upper = (bucket_num == procs - 1) ? INT_MAX : splitters[bucket_num];
    }

    // 6) All-to-all exchange of counts (small)

    delete[] int_samples;
    delete[] all_samples;
    delete[] splitters;
    delete[] buckets;
    break;
  }
  case DOUBLE:
    local_sort_mergesort(double_local_data, local_size, cfg.data_type);
    break;
  }
  CALI_MARK_END("comp");
  CALI_MARK_END("comp_large");

  // 2) Local sampling: pick s=P evenly spaced samples
  const int samples = procs;
  int *int_samples = nullptr;
  double *double_samples = nullptr;

  switch (cfg.data_type) {
  case INT:
    delete[] int_local_data;
    break;
  case DOUBLE:
    delete[] double_local_data;
    break;
  }

  delete[] dist.counts;
  delete[] dist.displs;

  return 0;
}

void merge_runs(const void *local_array, const void *tmp,
                const DATA_TYPE d_type, size_t left, size_t mid, size_t right) {
  int *int_local_array, int_tmp = nullptr;
  double *double_local_data, double_tmp = nullptr;

  size_t i = left, j = mid, k = left;
  switch (d_type) {
  case INT:
    int_local_array = (int *)local_array;
    int_tmp = (int *)tmp;
    while (i < mid && j < right)
      int_tmp[k++] = (int_local_array[i] <= int_local_array[j])
                         ? int_local_array[i++]
                         : int_local_array[j++];
    while (i < mid)
      int_tmp[k++] = int_local_array[i++];
    while (j < right)
      int_tmp[k++] = int_local_array[j++];
    for (size_t t = left; t < right; ++t)
      int_local_array[t] = int_tmp[t];
    break;
  case DOUBLE:
    double_local_data = (double *)local_array;
    double_tmp = (double *)tmp;
    while (i < mid && j < right)
      double_tmp[k++] = (double_local_array[i] <= double_local_array[j])
                            ? double_local_array[i++]
                            : double_local_array[j++];
    while (i < mid)
      double_tmp[k++] = double_local_array[i++];
    while (j < right)
      double_tmp[k++] = double_local_array[j++];
    for (size_t t = left; t < right; ++t)
      double_local_array[t] = double_tmp[t];
    break;
  }
}

void local_sort_mergesort(const void *local_array, const long long local_size,
                          const DATA_TYPE d_type) {
  int *int_tmp = nullptr;
  double *double_tmp = nullptr;

  switch (d_type) {
  case INT:
    int_tmp = new int[local_size];
    for (size_t width = 1; width < local_size; width <<= 1) {
      for (size_t left = 0; left < local_size; left += (width << 1)) {
        size_t mid = std::min(left + width, local_size);
        size_t right = std::min(left + (width << 1), local_size);
        if (mid < right)
          merge_runs(local_array, int_tmp, d_type, tmp, left, mid, right);
      }
    }
    delete[] int_tmp;
    break;
  case DOUBLE:
    double_tmp = new double[local_size];
    for (size_t width = 1; width < local_size; width <<= 1) {
      for (size_t left = 0; left < local_size; left += (width << 1)) {
        size_t mid = std::min(left + width, local_size);
        size_t right = std::min(left + (width << 1), local_size);
        if (mid < right)
          merge_runs(local_array, double_tmp, d_type, tmp, left, mid, right);
      }
    }
    delete[] double_tmp;
    break;
  }
}
