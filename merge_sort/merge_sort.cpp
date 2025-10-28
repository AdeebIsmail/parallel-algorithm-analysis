#include "mpi.h"
#include <adiak.hpp>
#include <caliper/cali-manager.h>
#include <caliper/cali.h>
#include <math.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

enum SortLevel { SORTED, PERTURBED, RANDOM, REVERSED };

// generate data in each process
// do comparison check after each merge sort,

template <typename T> T *merge_arrays(T *A, int nA, T *B, int nB) {
  T *result = new T[nA + nB];
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

template <typename T> int validate(T *local_data, int n) {
  for (int i = 1; i < n; i++) {
    if (local_data[i] < local_data[i - 1]) {
      return -1;
    }
  }
  return 1;
}

template <typename T> void MergeSort(T *arr, int n) {
  std::sort(arr, arr + n);
  // for (int current_size = 1; current_size < n;
  //      current_size = 2 * current_size) {
  //   for (int left_start = 0; left_start < n - 1;
  //        left_start += 2 * current_size) {
  //     int mid = left_start + current_size - 1;
  //     int right_end = left_start + 2 * current_size - 1;

  //     if (mid >= n)
  //       break;
  //     if (right_end >= n)
  //       right_end = n - 1;

  //     int left_size = mid - left_start + 1;
  //     int right_size = right_end - mid;

  //     int *left_arr = new int[left_size];
  //     int *right_arr = new int[right_size];

  //     for (int i = 0; i < left_size; i++)
  //       left_arr[i] = arr[left_start + i];
  //     for (int j = 0; j < right_size; j++)
  //       right_arr[j] = arr[mid + 1 + j];

  //     int i = 0, j = 0, k = left_start;

  //     while (i < left_size && j < right_size) {
  //       if (left_arr[i] <= right_arr[j]) {
  //         arr[k++] = left_arr[i++];
  //       } else {
  //         arr[k++] = right_arr[j++];
  //       }
  //     }

  //     while (i < left_size) {
  //       arr[k++] = left_arr[i++];
  //     }
  //     while (j < right_size) {
  //       arr[k++] = right_arr[j++];
  //     }

  //     delete[] left_arr;
  //     delete[] right_arr;
  //   }
  // }
}

template <typename T>
void run_merge_sort(int element, int sort_level, int rank, int size,
                    cali::ConfigManager &mgr, MPI_Datatype mpi_type,
                    const char *type_name) {
  int n;
  T *data = nullptr;
  n = 1 << element;
  int elements_per_proc = n / size;
  srand(101);

  CALI_MARK_BEGIN("data_init_runtime");
  data = new T[elements_per_proc];

  unsigned int seed = 12345 + rank * 1000;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> random_numbers(0, GLOBAL_MAX_NUMBER);
  std::uniform_int_distribution<> indexes(0, elements_per_proc - 1);

  if (sort_level == RANDOM) {
    for (int i = 0; i < elements_per_proc; i++) {
      data[i] = static_cast<T>(random_numbers(gen));
    }
  } else if (sort_level == SORTED) {
    int global_start = rank * elements_per_proc;
    int index = 0;
    for (int i = global_start; i < global_start + elements_per_proc; i++) {
      data[index] = static_cast<T>(i);
      index += 1;
    }
  } else if (sort_level == REVERSED) {
    int global_start = (size - rank - 1) * elements_per_proc;
    int index = 0;
    for (int i = global_start + elements_per_proc - 1; i >= global_start; i--) {
      data[index] = static_cast<T>(i);
      index += 1;
    }
  } else if (sort_level == PERTURBED) {
    int global_start = rank * elements_per_proc;
    int index = 0;
    for (int i = global_start; i < global_start + elements_per_proc; i++) {
      data[index] = static_cast<T>(i);
      index += 1;
    }
    int num_swaps = elements_per_proc / 100;
    for (int s = 0; s < num_swaps; s++) {
      int idx = rand() % (elements_per_proc - 1);
      int temp = data[idx];
      data[idx] = data[idx + 1];
      data[idx + 1] = temp;
    }
  }

  CALI_MARK_END("data_init_runtime");

  int local_n = n / size;
  T *local_data = data;

  CALI_MARK_BEGIN("comp");
  CALI_MARK_BEGIN("comp_small");
  MergeSort(local_data, local_n);
  CALI_MARK_END("comp_small");
  CALI_MARK_END("comp");

  int step = 1;
  while (step < size) {
    if (rank % (2 * step) == 0) {
      int recv_proc = rank + step;
      int recv_n;
      CALI_MARK_BEGIN("comm");
      CALI_MARK_BEGIN("comm_large");
      MPI_Recv(&recv_n, 1, MPI_INT, recv_proc, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      T *recv_data = new T[recv_n];
      MPI_Recv(recv_data, recv_n, mpi_type, recv_proc, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      CALI_MARK_END("comm_large");
      CALI_MARK_END("comm");
      CALI_MARK_BEGIN("comp");
      CALI_MARK_BEGIN("comp_large");
      T *merged = merge_arrays(local_data, local_n, recv_data, recv_n);
      CALI_MARK_END("comp_large");
      CALI_MARK_END("comp");
      CALI_MARK_BEGIN("correctness_check");
      if (validate(merged, local_n + recv_n) == -1) {
        cout << "Validation Failed" << endl;
        assert(false);
        delete[] local_data;
        delete[] recv_data;
        delete[] merged;
        return;
      }
      CALI_MARK_END("correctness_check");
      delete[] local_data;
      delete[] recv_data;
      local_data = merged;
      local_n += recv_n;

    } else {
      int send_proc = rank - step;
      CALI_MARK_BEGIN("comm");
      CALI_MARK_BEGIN("comm_large");
      MPI_Send(&local_n, 1, MPI_INT, send_proc, 0, MPI_COMM_WORLD);
      MPI_Send(local_data, local_n, mpi_type, send_proc, 0, MPI_COMM_WORLD);
      CALI_MARK_END("comm_large");
      CALI_MARK_END("comm");
      delete[] local_data;
      return;
    }
    step *= 2;
  }
  CALI_MARK_BEGIN("correctness_check");
  if (validate(local_data, local_n) == -1) {
    cout << "Validation Failed" << endl;
    assert(false);
  }
  CALI_MARK_END("correctness_check");

  delete[] local_data;

  adiak::init(NULL);
  adiak::launchdate();
  adiak::libraries();
  adiak::cmdline();
  adiak::clustername();
  adiak::value("algorithm", "Merge Sort");
  adiak::value("programming_model", "MPI");
  adiak::value("data_type", type_name);
  adiak::value("size_of_data_type", sizeof(T));
  adiak::value("input_size", n);
  adiak::value("input_type", sort_level);
  adiak::value("num_procs", size);
  adiak::value("group_num", "8");
  adiak::value("implementation_source",
               "https://stanford.edu/~rezab/classes/cme323/S16/notes/Lecture03/"
               "cme323_lec3.pdf");

  mgr.stop();
  mgr.flush();
}

int main(int argc, char **argv) {
  CALI_CXX_MARK_FUNCTION;

  int rank, size;
  int element;
  int sort_level;
  int use_float = 0; // 0 = int, 1 = float

  if (argc == 3) {
    element = atoi(argv[1]);
    sort_level = atoi(argv[2]);
  } else if (argc == 4) {
    element = atoi(argv[1]);
    sort_level = atoi(argv[2]);
    use_float = atoi(argv[3]);
  } else {
    printf(
        "\n Missing parameters. Usage: %s <element> <sort_level> [use_float]\n",
        argv[0]);
    printf(" element: power of 2 for array size\n");
    printf(" sort_level: 0=SORTED, 1=PERTURBED, 2=RANDOM, 3=REVERSED\n");
    printf(" use_float (optional): 0=int (default), 1=float\n");
    return 0;
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  cali::ConfigManager mgr;
  mgr.start();

  if (use_float == 1) {
    run_merge_sort<float>(element, sort_level, rank, size, mgr, MPI_FLOAT,
                          "float");
  } else {
    run_merge_sort<int>(element, sort_level, rank, size, mgr, MPI_INT, "int");
  }

  MPI_Finalize();
  return 0;
}
