#include <caliper/cali-manager.h>
#include <caliper/cali.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
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

  if (argc == 2) {
    element = atoi(argv[1]);
  } else {
    printf("\n Please provide how many elements");
    return 0;
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int n;
  int *data = nullptr;

  if (rank == 0) {
    CALI_MARK_BEGIN(data_init_runtime);
    n = 1 << element;
    data = new int[n];
    srand(101);
    for (int i = 0; i < n; i++) {
      data[i] = rand() % 100;
    }
    CALI_MARK_END(data_init_runtime);
  }
  // Size of the array
  CALI_MARK_BEGIN("comm");
  CALI_MARK_BEGIN("comm_small_Bcast");
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  CALI_MARK_END("comm_small_Bcast");
  CALI_MARK_END("comm");

  int local_n = n / size;
  int *local_data = new int[local_n];

  // Split up data of size local_n and send to each process, receives into
  // local_data and is of size local_n
  CALI_MARK_BEGIN("comm");
  CALI_MARK_BEGIN("comm_large_Scatter");
  MPI_Scatter(data, local_n, MPI_INT, local_data, local_n, MPI_INT, 0,
              MPI_COMM_WORLD);
  CALI_MARK_END("comm_large_Scatter");
  CALI_MARK_END("comm");

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

  if (rank == 0) {
    CALI_MARK_BEGIN(correctness_check);

    CALI_MARK_END(correctness_check);
  }

  delete[] local_data;
  if (rank == 0)
    delete[] data;

  MPI_Finalize();
  return 0;
}
