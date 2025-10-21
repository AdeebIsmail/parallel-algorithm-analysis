# CSCE 435 Group project

## 0. Group number and the communication medium(slack etc.):

Group number : 6
Communication medium : Discord

## 1. Group members, what algorithm they are implementing

1. Bitonic Sort: Jack Williams
2. Merge Sort: Adeeb Ismail
3. Radix Sort: Soham Nagawanshi
4. Sample Sort: David Armijo

### 2. Brief project description, what architecture you are comparing your sorting algorithms on.

### 2a. Pseudocode for Bitonic Sort.

BITONIC_MPI_SORT(comm, local[]):
    MPI_Init()
    P ← MPI_Comm_size(comm)
    r ← MPI_Comm_rank(comm)

    # 0) Local pre-sort (ascending)
    sort(local)   # e.g., quicksort; in-place; ascending

    # 1) Process-level bitonic network
    size ← 2
    while size ≤ P:

        # j walks partner distances inside this size-group: size/2, size/4, ..., 1
        j ← size / 2
        while j ≥ 1:

            partner ← r XOR j

            # Keep-lower vs keep-upper decision:
            # If our 'size' bit in r is 0 ⇒ we are in the group's lower half ⇒ keep lower n.
            # If it's 1 ⇒ we are in the group's upper half ⇒ keep upper n.
            keep_lower ← ((r AND size) == 0)

            # Exchange equal-sized sorted blocks with partner
            temp[] ← buffer of length n (same type as local[])
            MPI_Sendrecv(sendbuf=local, sendcount=n, sendtype=T, dest=partner, tag=0,
                         recvbuf=temp,  recvcount=n, recvtype=T, source=partner, recvtag=0,
                         comm, status=IGNORE)

            # Compare–split: merge two sorted blocks and keep our half
            if keep_lower:
                local ← MERGE_KEEP_LOW_N(local, temp)     # keep smallest n
            else:
                local ← MERGE_KEEP_HIGH_N(local, temp)    # keep largest n

            j ← j / 2
        end while

        size ← size * 2
    end while

    # (Optional) Gather final globally sorted data to root for verification/display
    # MPI_Gatherv(local → root)

    MPI_Finalize()
    return

# ----------------------- Helpers -----------------------

MERGE_KEEP_LOW_N(A[], B[]):
    # A and B: each length n, both sorted ascending
    # Return the smallest n elements of their union without building 2n
    i ← 0; j ← 0; k ← 0
    C[] ← new array length n
    while k < n:
        if j == n or (i < n and A[i] ≤ B[j]):
            C[k] ← A[i];  i ← i + 1
        else:
            C[k] ← B[j];  j ← j + 1
        k ← k + 1
    return C

MERGE_KEEP_HIGH_N(A[], B[]):
    # A and B: each length n, both sorted ascending
    # Return the largest n elements; scan from the end
    i ← n - 1; j ← n - 1; k ← n - 1
    C[] ← new array length n
    while k ≥ 0:
        if j < 0 or (i ≥ 0 and A[i] ≥ B[j]):
            C[k] ← A[i];  i ← i - 1
        else:
            C[k] ← B[j];  j ← j - 1
        k ← k - 1
    return C


Include MPI calls you will use to coordinate between processes.

### 2b. Pseudocode for Merge Sort.

Include MPI calls you will use to coordinate between processes.

#### The MPI calls we will need are

- MPI_Send
- MPI_Recv
- MPI_Init
- MPI_Finalize

#### Pseudocode (Master/Worker, Point-to-Point)

    IF master:
      for worker in workers:
        MPI_Send() send equal amounts of data to each process
      MergeSort() array of master
      for worker in workers:
        MPI_Recv() from all the workers
        Merge() the recieved arrays together
    ELSE:
      MPI_Recv() recieve data from the master
      MergeSort() array of each worker
      MPI_Send() send the sorted array to the master


    Merge(a,b) merges the two arrays together in sorted order
    MergeSort() performs merge sort

### 2c. Pseudocode for Radix Sort.

note: LSD = least significant digit
```
data: data to sort
k: maximum number of digits in any element of data
P: processor array

function radix-sort(data, k, P)
  # initialization
  blocks <- divide into |P| continguous blocks (each has size |data|/|P|)

  If Master: 
    MPI_Send(blocks_i, p_i) # send each block to corresponding processor
  
  If Worker:
    MPI_Recv(my_block, Master)
    Repeat for i in {1 ... k}
      buckets <- {0: [], 1: [], ..., 9: []} # itemizes data by i'th LSD
      hist <- [0, 0, .., 0] # keeps track of the size of each digit bucket
      for each element in my_block
        dig <- ith LSD of element
        add element to buckets[dig]
        hist[dig] += 1
      MPI_Allgather(hist) # broadcast histogram to all workers, result is |P| x 10 matrix
      # send buckets to processors in sequential fashion (i.e. bucket 0 will go to smaller rank processors but
                                                          bucket 9 will go to bigger rank processor) 
      recv_procs <- calculate which processes you will receive from based on hist
      send_procs <- calculate which processes you will send to based on hist
      for each recv_proc:
          MPI_Recv(blocks, recv_procs)
      for each send_proc:
          MPI_Send(blocks, send_procs)
    MPI_Send(blocks, Master)
  If Master:
    sorted_data <- empty array
    for each p in P: # iterate sequentially for correct order
        MPI_Recv(blocks, p)
        place blocks in sorted_data
      
    output sorted data
```
  


### 2d. Pseudocode for Sample Sort.

#### MPI Calls We will Need
- MPI_Send (For Master and Slave). Whenever point-to-point communication is necessary
- MPI_Recv (For Master and Slave). Whenever point-to-point communication is necessary
- MPI_Init
- **MPI_Allgather**
- **MPI_Alltoall**
- **MPI_Alltoallv**
- MPI_Finalize

#### Pseudocode
Include MPI calls you will use to coordinate between processes.
For array size N and processor count P,
1. Pre-Sorting (Sequential portion):
  - Split up and distribute the array evenly to each processor.
  - P1 gets matrix [0:1*(N/P)-1], P2 gets matrix [1*(N/P):2*(N/P)-1], ..., Pk gets matrix [(k-1)*(N/P):k*(N/P)-1] and so on.
  - Use MPI_Send to send the separated data to all the processes
2. Local Sort:
  - In each of the processors in the rank, sort the arrays locally.
3. Local Sampling:
  - In each of the local arrays in the rank, select `s = P` evenly spaced elements as its sample.
  - Example: `s=P=4`, positions samples will be (0,2,4,6) in each of those local arrays if the local array size is 8 elements
4. Gathering Samples:
  - Using `MPI_Allgather` function, gather all samples from every processor.
5. Sort the Gathered Samples:
  - Sort the gathered samples array (that came from all the ranks)
6. Partition by Splitters:
  - In the sorted gathered array, split that array up into P parts. In other words, make P-1 splitters and make them the values of your bucket list.
  - Ex: For 16 samples (from those 4 ranks or Processors), pick indices 1*(16/4), 2\*(16/4), and 3\*(16/4). These will be your S1 = __, S2 = __, and S3 = __.
7. Partition by splitters.
  - Your buckets are the values at the indices of the sorted, combined sample list. All the ranks then receive that list of splitter values to bin their local lists via `MPI_Alltoall`.
8. Sort local lists using global splitter values.
  - Make sure to place them in their respective bins if they are less than or equal to that splitter value.
9. All-to-All Exchange
  - Each rank sends Bj to rank j. After `MPI_Alltoallv`, each rank holds one bucket's *global* range.
10. Final local sort/merge:
  - If you've done stable partitioning on already sorted local arrays, each rank can *merge* its received, already-sorted chunks (k-way merge). Otherwise, just sort once more locally. The concatenation across (0...3) is the globally sorted order

### 3. Evaluation plan - what and how will you measure and compare

We will vary multiple parameters and will measure performance
using caliper instrumentation. For each experiment, we will take the measurements
mentioned in 6c. Our first experiment will vary data types and fix the
other parameters. For this experiment, we plan to compare sorting an integer array with a double
array. Our second experiment will vary input data sizes while fixing other parameters. We will use
2^16, 2^18, 2^20, 2^22, 2^24, 2^26, and 2^28 elements for evaluation. Our third experiment will
vary the sorting level. We will test on a sorted array, a sorted array with 1% perturbation,
a randomly generated array, and a reversely sorted array. To gauge strong scaling performance,
we will fix the input data size (likely to 2^22 elements) and use
2, 4, 8, 16, 32, 64, 128, 256, 512, and 1024 processors for evaluation. To gauge weak scaling,
we will vary the input data size as we increase the processor count. We will likely use
the following pairs:

* (2^16 elements, 16 processors)
* (2^18 elements, 32 processors)
* (2^20 elements, 64 processors)
* (2^22 elements, 128 processors)
* (2^24 elements, 256 processors)
* (2^26 elements, 512 processors)
* (2^28 elements, 1024 processors)

We will carry out these experiments for each of the sorting algorithms
and compare the performance across different input sizes / processor counts using
graphics. We will also include graphics that compare performance between algorithms.

### 4. Caliper instrumentation

Please use the caliper build `/scratch/group/csce-435-f25/Caliper/caliper/share/cmake/caliper`
(same as lab2 build.sh) to collect caliper files for each experiment you run.

Your Caliper annotations should result in the following calltree
(use `Thicket.tree()` to see the calltree):

```
main
|_ data_init_X      # X = runtime OR io
|_ comm
|    |_ comm_small
|    |_ comm_large
|_ comp
|    |_ comp_small
|    |_ comp_large
|_ correctness_check
```

Required region annotations:

- `main` - top-level main function.
  - `data_init_X` - the function where input data is generated or read in from file. Use _data_init_runtime_ if you are generating the data during the program, and _data_init_io_ if you are reading the data from a file.
  - `correctness_check` - function for checking the correctness of the algorithm output (e.g., checking if the resulting data is sorted).
  - `comm` - All communication-related functions in your algorithm should be nested under the `comm` region.
    - Inside the `comm` region, you should create regions to indicate how much data you are communicating (i.e., `comm_small` if you are sending or broadcasting a few values, `comm_large` if you are sending all of your local values).
    - Notice that auxillary functions like MPI_init are not under here.
  - `comp` - All computation functions within your algorithm should be nested under the `comp` region.
    - Inside the `comp` region, you should create regions to indicate how much data you are computing on (i.e., `comp_small` if you are sorting a few values like the splitters, `comp_large` if you are sorting values in the array).
    - Notice that auxillary functions like data_init are not under here.
  - `MPI_X` - You will also see MPI regions in the calltree if using the appropriate MPI profiling configuration (see **Builds/**). Examples shown below.

All functions will be called from `main` and most will be grouped under either `comm` or `comp` regions, representing communication and computation, respectively. You should be timing as many significant functions in your code as possible. **Do not** time print statements or other insignificant operations that may skew the performance measurements.

### **Nesting Code Regions Example** - all computation code regions should be nested in the "comp" parent code region as following:

```
CALI_MARK_BEGIN("comp");
CALI_MARK_BEGIN("comp_small");
sort_pivots(pivot_arr);
CALI_MARK_END("comp_small");
CALI_MARK_END("comp");

# Other non-computation code
...

CALI_MARK_BEGIN("comp");
CALI_MARK_BEGIN("comp_large");
sort_values(arr);
CALI_MARK_END("comp_large");
CALI_MARK_END("comp");
```

### **Calltree Example**:

```
# MPI Mergesort
4.695 main
├─ 0.001 MPI_Comm_dup
├─ 0.000 MPI_Finalize
├─ 0.000 MPI_Finalized
├─ 0.000 MPI_Init
├─ 0.000 MPI_Initialized
├─ 2.599 comm
│  ├─ 2.572 MPI_Barrier
│  └─ 0.027 comm_large
│     ├─ 0.011 MPI_Gather
│     └─ 0.016 MPI_Scatter
├─ 0.910 comp
│  └─ 0.909 comp_large
├─ 0.201 data_init_runtime
└─ 0.440 correctness_check
```

### 5. Collect Metadata

Have the following code in your programs to collect metadata:

```
adiak::init(NULL);
adiak::launchdate();    // launch date of the job
adiak::libraries();     // Libraries used
adiak::cmdline();       // Command line used to launch the job
adiak::clustername();   // Name of the cluster
adiak::value("algorithm", algorithm); // The name of the algorithm you are using (e.g., "merge", "bitonic")
adiak::value("programming_model", programming_model); // e.g. "mpi"
adiak::value("data_type", data_type); // The datatype of input elements (e.g., double, int, float)
adiak::value("size_of_data_type", size_of_data_type); // sizeof(datatype) of input elements in bytes (e.g., 1, 2, 4)
adiak::value("input_size", input_size); // The number of elements in input dataset (1000)
adiak::value("input_type", input_type); // For sorting, this would be choices: ("Sorted", "ReverseSorted", "Random", "1_perc_perturbed")
adiak::value("num_procs", num_procs); // The number of processors (MPI ranks)
adiak::value("scalability", scalability); // The scalability of your algorithm. choices: ("strong", "weak")
adiak::value("group_num", group_number); // The number of your group (integer, e.g., 1, 10)
adiak::value("implementation_source", implementation_source); // Where you got the source code of your algorithm. choices: ("online", "ai", "handwritten").
```

They will show up in the `Thicket.metadata` if the caliper file is read into Thicket.

### **See the `Builds/` directory to find the correct Caliper configurations to get the performance metrics.** They will show up in the `Thicket.dataframe` when the Caliper file is read into Thicket.

## 6. Performance evaluation

Include detailed analysis of computation performance, communication performance.
Include figures and explanation of your analysis.

### 6a. Vary the following parameters

For input_size's:

- 2^16, 2^18, 2^20, 2^22, 2^24, 2^26, 2^28

For input_type's:

- Sorted, Random, Reverse sorted, 1%perturbed

MPI: num_procs:

- 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024

This should result in 4x7x10=280 Caliper files for your MPI experiments.

### 6b. Hints for performance analysis

To automate running a set of experiments, parameterize your program.

- input_type: "Sorted" could generate a sorted input to pass into your algorithms
- algorithm: You can have a switch statement that calls the different algorithms and sets the Adiak variables accordingly
- num_procs: How many MPI ranks you are using

When your program works with these parameters, you can write a shell script
that will run a for loop over the parameters above (e.g., on 64 processors,
perform runs that invoke algorithm2 for Sorted, ReverseSorted, and Random data).

### 6c. You should measure the following performance metrics

- `Time`
  - Min time/rank
  - Max time/rank
  - Avg time/rank
  - Total time
  - Variance time/rank

## 7. Presentation

Plots for the presentation should be as follows:

- For each implementation:
  - For each of comp_large, comm, and main:
    - Strong scaling plots for each input_size with lines for input_type (7 plots - 4 lines each)
    - Strong scaling speedup plot for each input_type (4 plots)
    - Weak scaling plots for each input_type (4 plots)

Analyze these plots and choose a subset to present and explain in your presentation.

## 8. Final Report

Submit a zip named `TeamX.zip` where `X` is your team number. The zip should contain the following files:

- Algorithms: Directory of source code of your algorithms.
- Data: All `.cali` files used to generate the plots seperated by algorithm/implementation.
- Jupyter notebook: The Jupyter notebook(s) used to generate the plots for the report.
- Report.md
