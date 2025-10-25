#include <adiak.hpp>
#include <algorithm>
#include <caliper/cali-manager.h>
#include <caliper/cali.h>
#include <cassert>
#include <math.h>
#include <mpi.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

enum SortLevel { SORTED, PERTURBED, RANDOM, REVERSED };
#define GLOBAL_MAX_NUMBER 400e6
// generate data in each process
// do comparison check after each merge sort,
int *merge_arrays(int *A, int nA, int *B, int nB) {
  int *result = new int[nA + nB];
  int i = 0, j = 0, k = 0;

  // Merge the two arrays
  while (i < nA && j < nB) {
    if (A[i] <= B[j]) {
      result[k++] = A[i++];
    } else {
      result[k++] = B[j++];
    }
  }

  while (i < nA) {
    result[k++] = A[i++];
  }

  while (j < nB) {
    result[k++] = B[j++];
  }

  return result;
}

int validate(int *local_data, int n) {
  for (int i = 1; i < n; i++) {
    if (local_data[i] < local_data[i - 1]) {
      return -1;
    }
  }
  return 1;
}

void MergeSort(int *arr, int n) {

  for (int current_size = 1; current_size < n;
       current_size = 2 * current_size) {
    for (int left_start = 0; left_start < n - 1;
         left_start += 2 * current_size) {
      int mid = left_start + current_size - 1;
      int right_end = left_start + 2 * current_size - 1;

      if (mid >= n)
        break;
      if (right_end >= n)
        right_end = n - 1;

      int left_size = mid - left_start + 1;
      int right_size = right_end - mid;

      int *left_arr = new int[left_size];
      int *right_arr = new int[right_size];

      for (int i = 0; i < left_size; i++)
        left_arr[i] = arr[left_start + i];
      for (int j = 0; j < right_size; j++)
        right_arr[j] = arr[mid + 1 + j];

      int i = 0, j = 0, k = left_start;

      while (i < left_size && j < right_size) {
        if (left_arr[i] <= right_arr[j]) {
          arr[k++] = left_arr[i++];
        } else {
          arr[k++] = right_arr[j++];
        }
      }

      while (i < left_size) {
        arr[k++] = left_arr[i++];
      }
      while (j < right_size) {
        arr[k++] = right_arr[j++];
      }

      delete[] left_arr;
      delete[] right_arr;
    }
  }
}

int main(int argc, char **argv) {
  CALI_CXX_MARK_FUNCTION;

  int rank, size;
  int element;
  int sort_level;

  if (argc == 3) {
    element = atoi(argv[1]);
    sort_level = atoi(argv[2]);

  } else {
    printf("\n Missing parameters");
    return 0;
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  cali::ConfigManager mgr;
  mgr.start();

  int n;
  int *data = nullptr;
  n = 1 << element;
  int elements_per_proc = n / size;

  CALI_MARK_BEGIN("data_init_runtime");
  data = new int[elements_per_proc];

  unsigned int seed = 12345 + rank * 1000;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> random_numbers(0, GLOBAL_MAX_NUMBER);
  std::uniform_int_distribution<> indexes(0, elements_per_proc - 1);

  if (sort_level == RANDOM) {
    for (int i = 0; i < elements_per_proc; i++) {
      data[i] = random_numbers(gen);
    }
  } else if (sort_level == SORTED) {
    int global_start = rank * elements_per_proc;
    int index = 0;
    for (int i = global_start; i < global_start + elements_per_proc; i++) {
      data[index] = i;
      index += 1;
    }
  } else if (sort_level == REVERSED) {
    int global_start = (size - rank - 1) * elements_per_proc;
    int index = 0;
    for (int i = global_start + elements_per_proc - 1; i >= global_start; i--) {
      data[index] = i;
      index += 1;
    }
  } else if (sort_level == PERTURBED) {
    int global_start = rank * elements_per_proc;
    int index = 0;
    for (int i = global_start; i < global_start + elements_per_proc; i++) {
      data[index] = i;
      index += 1;
    }
    int num_swaps = elements_per_proc / 100;
    for (int s = 0; s < num_swaps; s++) {
      int idx1 = indexes(gen);
      int idx2 = indexes(gen);
      while (idx1 == idx2 && elements_per_proc > 1) {
        idx2 = indexes(gen);
      }
      std::swap(data[idx1], data[idx2]);
    }
  }

  CALI_MARK_END("data_init_runtime");

  int local_n = n / size;
  int *local_data = data;

  for (int i = 0; i < local_n; i++) {
    cout << local_data[i] << " ";
  }
  cout << "" << endl;

  CALI_MARK_BEGIN("comp");
  CALI_MARK_BEGIN("comp_small_merge_sort");
  MergeSort(local_data, local_n);
  CALI_MARK_END("comp_small_merge_sort");
  CALI_MARK_END("comp");

  int step = 1;
  while (step < size) {
    if (rank % (2 * step) == 0) {
      int recv_proc = rank + step;
      int recv_n;
      CALI_MARK_BEGIN("comm");
      CALI_MARK_BEGIN("comm_large_recv");
      MPI_Recv(&recv_n, 1, MPI_INT, recv_proc, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      int *recv_data = new int[recv_n];
      MPI_Recv(recv_data, recv_n, MPI_INT, recv_proc, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      CALI_MARK_END("comm_large_recv");
      CALI_MARK_END("comm");
      CALI_MARK_BEGIN("comp");
      CALI_MARK_BEGIN("comp_large_merge_arrays");
      int *merged = merge_arrays(local_data, local_n, recv_data, recv_n);
      CALI_MARK_END("comp_large_merge_arrays");
      CALI_MARK_END("comp");
      CALI_MARK_BEGIN("correctness_check");
      if (validate(merged, local_n + recv_n) == -1) {
        cout << "Validation Failed" << endl;
        assert(false);
        break;
      }
      CALI_MARK_END("correctness_check");
      delete[] local_data;
      delete[] recv_data;
      local_data = merged;
      local_n += recv_n;

    } else {
      int send_proc = rank - step;
      CALI_MARK_BEGIN("comm");
      CALI_MARK_BEGIN("comm_large_send");
      MPI_Send(&local_n, 1, MPI_INT, send_proc, 0, MPI_COMM_WORLD);
      MPI_Send(local_data, local_n, MPI_INT, send_proc, 0, MPI_COMM_WORLD);
      CALI_MARK_END("comm_large_send");
      CALI_MARK_END("comm");
      break;
    }
    step *= 2;
  }
  cout << "ENding" << endl;
  if (rank == 0) {
    for (int i = 0; i < local_n; i++) {
      cout << local_data[i] << " ";
    }
  }

  delete[] local_data;

  MPI_Finalize();

  adiak::init(NULL);
  adiak::launchdate();  // launch date of the job
  adiak::libraries();   // Libraries used
  adiak::cmdline();     // Command line used to launch the job
  adiak::clustername(); // Name of the cluster
  adiak::value("algorithm", "Merge Sort");  // The name of the algorithm you are
                                            // using (e.g., "merge", "bitonic")
  adiak::value("programming_model", "MPI"); // e.g. "mpi"
  adiak::value(
      "data_type",
      "int"); // The datatype of input elements (e.g., double, int, float)
  adiak::value("size_of_data_type",
               sizeof(int)); // sizeof(datatype) of input elements in
                             // bytes (e.g., 1, 2, 4)
  adiak::value("input_size",
               n); // The number of elements in input dataset (1000)
  adiak::value("input_type",
               sort_level); // For sorting, this would be choices: ("Sorted",
                            // "ReverseSorted", "Random", "1_perc_perturbed")
  adiak::value("num_procs", size); // The number of processors (MPI ranks)
  // adiak::value("scalability",
  //              "scalability"); // The scalability of your algorithm. choices:
  //                            // ("strong", "weak")
  adiak::value("group_num",
               "8"); // The number of your group (integer, e.g., 1, 10)
  adiak::value("implementation_source",
               "https://stanford.edu/~rezab/classes/cme323/S16/notes/Lecture03/"
               "cme323_lec3.pdf"); // Where you got the source code of
                                   // your algorithm. choices:
                                   // ("online", "ai", "handwritten").

  mgr.stop();
  mgr.flush();

  return 0;
}
