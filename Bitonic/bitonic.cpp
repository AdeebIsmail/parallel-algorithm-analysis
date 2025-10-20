// bitonic_mpi.cpp  (per-rank generation + debug prints + fixed bitonic logic)
// Build (Intel MPI example):
//   mpiicpc -O3 bitonic_mpi.cpp -o bitonic_mpi -lcaliper -ladiak
// or OpenMPI:
//   mpicxx  -O3 bitonic_mpi.cpp -o bitonic_mpi -lcaliper -ladiak
// (Make sure C++17 is enabled if you use CMake: set(CMAKE_CXX_STANDARD 17))
// Run:
//   mpirun -np 4 ./bitonic_mpi --n 1024 --input_type Random --verify

#include <mpi.h>
#include <vector>
#include <algorithm>
#include <random>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cassert>

#include <caliper/cali.h>
#include <caliper/cali-manager.h>
#include <adiak.hpp>

// --------------------- Simple CLI ---------------------
struct Args {
    long long N = 1LL << 20;            // default = 2^20 elements
    std::string input_type = "Random";  // Sorted | ReverseSorted | Random | 1_perc_perturbed
    bool verify = false;
};
Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--n") && i+1 < argc)        a.N = std::stoll(argv[++i]);
        else if (!std::strcmp(argv[i], "--input_type") && i+1 < argc) a.input_type = argv[++i];
        else if (!std::strcmp(argv[i], "--verify"))            a.verify = true;
    }
    return a;
}

// ----------- Debug helpers -----------
static inline void rank_print_ordered(MPI_Comm comm, int r, int P, const std::string& s) {
    for (int t = 0; t < P; ++t) {
        MPI_Barrier(comm);
        if (r == t) { std::cout << s << std::flush; }
    }
    MPI_Barrier(comm);
}
static inline std::pair<int,int> minmax_of(const std::vector<int>& v) {
    return { v.front(), v.back() }; // assumes non-empty and sorted
}

// --------------------- Per-rank local data generation (with prints) ---------------------
void make_local_data(std::vector<int>& local, long long N, int r, int P, const Args& a) {
    CALI_CXX_MARK_SCOPE("data_init_runtime");

    const long long n = (long long)local.size();
    const long long start = (long long)r * n;  // global start index of this rank's slice

    if (a.input_type == "Sorted") {
        for (long long i = 0; i < n; ++i)
            local[(size_t)i] = (int)(start + i);

    } else if (a.input_type == "ReverseSorted") {
        for (long long i = 0; i < n; ++i) {
            long long g = start + i;                 // global index
            local[(size_t)i] = (int)(N - 1 - g);     // reverse of global index
        }

    } else if (a.input_type == "1_perc_perturbed") {
        // Start from the globally-consistent "Sorted" slice, then perturb ~1% within this slice.
        for (long long i = 0; i < n; ++i)
            local[(size_t)i] = (int)(start + i);
        std::mt19937 rng(42u + (unsigned)r); // per-rank seed (reproducible)
        std::uniform_int_distribution<long long> pos(0, n - 1);
        const long long swaps = std::max(1LL, n/100);
        for (long long k = 0; k < swaps; ++k) {
            long long i = pos(rng), j = pos(rng);
            std::swap(local[(size_t)i], local[(size_t)j]);
        }

    } else { // Random
        std::mt19937 rng(42u + (unsigned)r); // per-rank seed to avoid identical blocks
        std::uniform_int_distribution<int> dist(0, 1000000000);
        for (long long i = 0; i < n; ++i)
            local[(size_t)i] = dist(rng);
    }

    // ----- PRINT: what this rank generated -----
    {
        std::ostringstream oss;
        oss << "[gen] rank " << r << "/" << P
            << " N=" << N << " n=" << n
            << " type=" << a.input_type;

        if (!local.empty()) {
            auto minmax = std::minmax_element(local.begin(), local.end());
            oss << " local[min,max]=[" << *minmax.first << "," << *minmax.second << "]";
        }
        oss << "\n";

        // For tiny blocks, dump contents
        if ((int)n <= 32) {
            oss << "       rank " << r << " data: ";
            for (int i = 0; i < (int)n; ++i) oss << local[(size_t)i] << (i+1==(int)n?'\n':' ');
        }
        rank_print_ordered(MPI_COMM_WORLD, r, P, oss.str());
    }
}

// --------------------- Correctness check (with prints) ---------------------
// (1) Each local block ascending
// (2) Boundaries: max[r] <= min[r+1] for all r
bool is_sorted_global(const std::vector<int>& local, MPI_Comm comm) {
    CALI_CXX_MARK_SCOPE("correctness_check");

    int P, r;
    MPI_Comm_size(comm, &P);
    MPI_Comm_rank(comm, &r);

    // Local ascending check
    bool ok_local = std::is_sorted(local.begin(), local.end());
    int my_min = local.empty() ? 0 : local.front();
    int my_max = local.empty() ? 0 : local.back();

    // Local status print
    {
        std::ostringstream oss;
        oss << "[verify/local] rank " << r << "/" << P
            << " sorted=" << (ok_local ? "Y" : "N")
            << " min=" << my_min << " max=" << my_max << "\n";
        rank_print_ordered(comm, r, P, oss.str());
    }

    // Gather mins/maxs to rank 0 for boundary prints
    std::vector<int> mins, maxs;
    if (r == 0) { mins.resize(P); maxs.resize(P); }
    MPI_Gather(&my_min, 1, MPI_INT, mins.data(), 1, MPI_INT, 0, comm);
    MPI_Gather(&my_max, 1, MPI_INT, maxs.data(), 1, MPI_INT, 0, comm);

    int boundary_ok_all = 1;
    if (r == 0) {
        std::ostringstream oss;
        oss << "[verify/boundaries] checks (max[r] <= min[r+1]):\n";
        for (int i = 0; i < P - 1; ++i) {
            bool ok = (maxs[i] <= mins[i+1]);
            boundary_ok_all &= (ok ? 1 : 0);
            oss << "  r=" << i << " → r=" << (i+1)
                << " : " << maxs[i] << " <= " << mins[i+1]
                << "  => " << (ok ? "OK" : "FAIL") << "\n";
        }
        std::cout << oss.str() << std::flush;
    }

    int ok_local_all = ok_local ? 1 : 0, ok_local_reduce = 0;
    MPI_Allreduce(&ok_local_all, &ok_local_reduce, 1, MPI_INT, MPI_MIN, comm);
    MPI_Bcast(&boundary_ok_all, 1, MPI_INT, 0, comm);

    return (ok_local_reduce == 1) && (boundary_ok_all == 1);
}

// --------------------- Merge-keep helpers ---------------------
static inline void merge_keep_low_n(const std::vector<int>& A,
                                    const std::vector<int>& B,
                                    std::vector<int>& out) {
    CALI_CXX_MARK_FUNCTION;
    size_t n = out.size();
    size_t i=0, j=0, k=0;
    while (k < n) {
        if (j == n || (i < n && A[i] <= B[j])) out[k++] = A[i++];
        else                                   out[k++] = B[j++];
    }
}
static inline void merge_keep_high_n(const std::vector<int>& A,
                                     const std::vector<int>& B,
                                     std::vector<int>& out) {
    CALI_CXX_MARK_FUNCTION;
    size_t n = out.size();
    long long i = (long long)n - 1, j = (long long)n - 1, k = (long long)n - 1;
    while (k >= 0) {
        if (j < 0 || (i >= 0 && A[(size_t)i] >= B[(size_t)j])) out[(size_t)k--] = A[(size_t)i--];
        else                                                    out[(size_t)k--] = B[(size_t)j--];
    }
}

// --------------------- Bitonic process-level network ---------------------
void bitonic_sort_mpi(std::vector<int>& local, MPI_Comm comm) {
    //CALI_CXX_MARK_FUNCTION;

    int P, r;
    MPI_Comm_size(comm, &P);
    MPI_Comm_rank(comm, &r);

    const int n = (int)local.size();
    std::vector<int> tmp(n), next(n); // partner buf & output buf

    // 0) Local pre-sort
    {
        CALI_MARK_BEGIN("comp");
        CALI_MARK_BEGIN("comp_large");
        std::sort(local.begin(), local.end());
        CALI_MARK_END("comp_large");
        CALI_MARK_END("comp");
    }

    // 1) Process-level bitonic network
    for (int size = 2; size <= P; size <<= 1) {
        for (int j = size >> 1; j >= 1; j >>= 1) {
            const int partner = r ^ j;

            // Group direction: lower half ASC, upper half DESC
            const bool ascending_group = ((r & size) == 0);

            // Pair decision: when ascending, lower rank keeps LOWER half; when descending, higher keeps LOWER half
            const bool i_am_lower_rank = (r < partner);
            const bool keep_lower = ascending_group ? i_am_lower_rank
                                                    : !i_am_lower_rank;

            // BEFORE snapshot
            auto [pre_min, pre_max] = minmax_of(local);
            {
                std::ostringstream oss;
                oss << "[rank " << r << "] size=" << size
                    << " j=" << j
                    << " partner=" << partner
                    << " decision=" << (keep_lower ? "KEEP_LOW" : "KEEP_HIGH")
                    << " pre[min,max]=[" << pre_min << "," << pre_max << "]\n";
                rank_print_ordered(comm, r, P, oss.str());
            }

            // Exchange
            {
                CALI_MARK_BEGIN("comm");
                CALI_MARK_BEGIN("comm_large");
                MPI_Sendrecv(local.data(), n, MPI_INT, partner, 0,
                             tmp.data(),   n, MPI_INT, partner, 0,
                             comm, MPI_STATUS_IGNORE);
                CALI_MARK_END("comm_large");
                CALI_MARK_END("comm");
            }

            // Compare-split
            {
                CALI_MARK_BEGIN("comp");
                CALI_MARK_BEGIN("comp_large");
                if (keep_lower) merge_keep_low_n(local, tmp, next);
                else            merge_keep_high_n(local, tmp, next);
                local.swap(next);
                CALI_MARK_END("comp_large");
                CALI_MARK_END("comp");
            }

            // AFTER snapshot
            bool ok_local_after = std::is_sorted(local.begin(), local.end());
            auto [post_min, post_max] = minmax_of(local);
            {
                std::ostringstream oss;
                oss << "  [rank " << r << "] POST size=" << size
                    << " j=" << j
                    << " post[min,max]=[" << post_min << "," << post_max << "]"
                    << " sorted=" << (ok_local_after ? "Y" : "N") << "\n";
                rank_print_ordered(comm, r, P, oss.str());
            }
            if (!ok_local_after) {
                if (r == 0) std::cerr << "Local block not sorted after step. Aborting.\n";
                MPI_Abort(comm, 2);
            }
        }
    }
}

// --------------------- Main ---------------------
int main(int argc, char** argv) {
    CALI_CXX_MARK_FUNCTION;
    MPI_Init(&argc, &argv);

    // Adiak metadata init
    adiak::init(nullptr);
    adiak::launchdate();
    adiak::libraries();
    adiak::cmdline();
    adiak::clustername();

    // MPI world
    MPI_Comm comm = MPI_COMM_WORLD;
    int P = 0, r = 0;
    MPI_Comm_size(comm, &P);
    MPI_Comm_rank(comm, &r);

    // Args
    Args args = parse_args(argc, argv);

    // Problem setup
    const long long N = args.N;
    if (N % P != 0) {
        if (r == 0) std::cerr << "[Error] N must be divisible by P (got N=" << N << ", P=" << P << ").\n";
        MPI_Abort(comm, 1);
    }
    const int n = (int)(N / P);

    // Per-rank data generation (no root scatter)
    std::vector<int> local(n);
    make_local_data(local, N, r, P, args);

    // Adiak run metadata
    {
        std::string algorithm = "bitonic";
        std::string programming_model = "mpi";
        std::string data_type = "int";
        int size_of_data_type = sizeof(int);
        std::string input_type = args.input_type;
        int input_size = (int)N;
        int num_procs = P;
        std::string scalability = "strong";
        int group_number = 6;
        std::string implementation_source = "handwritten";

        adiak::value("algorithm", algorithm);
        adiak::value("programming_model", programming_model);
        adiak::value("data_type", data_type);
        adiak::value("size_of_data_type", size_of_data_type);
        adiak::value("input_size", input_size);
        adiak::value("input_type", input_type);
        adiak::value("num_procs", num_procs);
        adiak::value("scalability", scalability);
        adiak::value("group_num", group_number);
        adiak::value("implementation_source", implementation_source);
    }

    // Sort
    bitonic_sort_mpi(local, comm);

    // Verify & report
    bool ok = is_sorted_global(local, comm);
    if (r == 0) {
        std::cout << "[RESULT] bitonic_mpi verify=" << (ok ? "PASS" : "FAIL")
                  << " N=" << args.N << " P=" << P
                  << " input_type=" << args.input_type << "\n";
    }
    if (!ok) MPI_Abort(comm, 1);

    MPI_Finalize();
    return 0;
}
