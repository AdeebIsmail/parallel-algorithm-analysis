// main.cpp — MPI Sample Sort with Caliper/Adiak instrumentation (C++14-compatible)
// Build: mpicxx -O3 -std=c++14 main.cpp -lcaliper -ladiak
// Run:   mpirun -np 4 ./a.out --n 1048576 --datatype int --input random
//
// Notes:
// - No `if constexpr` used; code compiles under C++14+
// - Caliper regions: main, data_init_runtime, comp{small,large}, comm{small,large}, correctness_check

#include <mpi.h>
#include <caliper/cali-manager.h>
#include <caliper/cali.h>
#include <adiak.hpp>

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

// ----------------------------- CLI -----------------------------

struct Options {
  long long N = 1<<20;          // global element count
  std::string datatype = "int"; // "int" or "double"
  std::string input = "random"; // "random", "sorted", "reversed", "perturbed1"
  uint64_t seed = 42;
  int group_num = 6;
  std::string scalability = "strong";
};

static void usage(const char* prog) {
  if (!prog) prog = "sample_sort";
  std::cerr
    << "Usage: " << prog << " [--n N] [--datatype int|double] [--input random|sorted|reversed|perturbed1] [--seed S]\n"
    << "Examples:\n"
    << "  mpirun -np 8 " << prog << " --n 1048576 --datatype int --input random\n";
}

// Very small parser
static Options parse_cli(int argc, char** argv, int rank) {
  Options opt;
  for (int i=1;i<argc;++i) {
    std::string a = argv[i];
    auto need = [&](int i){ if (i+1>=argc){ if(!rank) usage(argv[0]); MPI_Abort(MPI_COMM_WORLD,1);} };
    if (a=="--n"){ need(i); opt.N = std::atoll(argv[++i]); }
    else if (a=="--datatype"){ need(i); opt.datatype = argv[++i]; }
    else if (a=="--input"){ need(i); opt.input = argv[++i]; }
    else if (a=="--seed"){ need(i); opt.seed = std::strtoull(argv[++i],nullptr,10); }
    else if (a=="--help" || a=="-h"){ if(!rank) usage(argv[0]); MPI_Abort(MPI_COMM_WORLD,0); }
    else if (a=="--group"){ need(i); opt.group_num = std::atoi(argv[++i]); }
    else if (a=="--scalability"){ need(i); opt.scalability = argv[++i]; }
    else {
      if (!rank) std::cerr << "Unknown arg: " << a << "\n";
      if (!rank) usage(argv[0]);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  return opt;
}

// ----------------------------- Utils -----------------------------

template<class T>
static bool is_sorted_non_decreasing(const std::vector<T>& v) {
  for (size_t i=1;i<v.size();++i)
    if (v[i] < v[i-1]) return false;
  return true;
}

// Correctness check without gathering all data.
// 1) Ensure each rank's local data is sorted.
// 2) Ensure boundary monotonicity: last element of rank r-1 <= first element of rank r.
template<class T>
static bool global_sorted_check(const std::vector<T>& local, MPI_Comm comm) {
  int rank, P; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&P);

  bool loc_sorted = is_sorted_non_decreasing(local);
  int loc_ok = loc_sorted ? 1 : 0;
  int all_ok = 0;
  MPI_Allreduce(&loc_ok, &all_ok, 1, MPI_INT, MPI_MIN, comm);
  if (!all_ok) return false;

  T first = (local.empty() ? std::numeric_limits<T>::lowest() : local.front());
  T last  = (local.empty() ? std::numeric_limits<T>::lowest() : local.back());

  // Exchange neighbor boundaries (two-way)
  MPI_Datatype Ttype = std::is_same<T,double>::value ? MPI_DOUBLE : MPI_LONG_LONG;

  T recv_prev_last = std::numeric_limits<T>::lowest();
  if (rank>0) {
    T send_first = first;
    MPI_Sendrecv(&send_first, 1, Ttype, rank-1, 11,
                 &recv_prev_last, 1, Ttype, rank-1, 12, comm, MPI_STATUS_IGNORE);
  }
  T recv_next_first = std::numeric_limits<T>::max();
  if (rank+1 < P) {
    T send_last = last;
    MPI_Sendrecv(&send_last, 1, Ttype, rank+1, 12,
                 &recv_next_first, 1, Ttype, rank+1, 11, comm, MPI_STATUS_IGNORE);
  }

  int boundary_ok = 1;
  if (rank>0 && !local.empty()) {
    if (recv_prev_last > first) boundary_ok = 0;
  }
  if (rank+1<P && !local.empty()) {
    if (last > recv_next_first && recv_next_first!=std::numeric_limits<T>::max()) boundary_ok = 0;
  }
  int gboundary_ok = 0;
  MPI_Allreduce(&boundary_ok, &gboundary_ok, 1, MPI_INT, MPI_MIN, comm);
  return gboundary_ok==1;
}

// K-way merge of sorted runs concatenated in 'in' with run sizes in 'sizes' (sum = total).
// Produces sorted 'out'.
template<class T>
static void kway_merge(const std::vector<T>& in,
                       const std::vector<int>& sizes,
                       std::vector<T>& out) {
  int K = (int)sizes.size();
  struct Run { size_t begin, end, idx; };
  std::vector<Run> runs; runs.reserve(K);
  size_t pos = 0;
  for (int k=0;k<K;++k) {
    runs.push_back({pos, pos + (size_t)sizes[k], pos});
    pos += (size_t)sizes[k];
  }
  out.clear(); out.reserve(in.size());
  struct Node { T v; int r; };
  struct Cmp { bool operator()(const Node& a, const Node& b) const { return a.v > b.v; } };
  std::vector<Node> heap;
  heap.reserve(K);
  for (int r=0;r<K;++r) if (runs[r].idx < runs[r].end) {
    heap.push_back({in[runs[r].idx], r});
  }
  std::make_heap(heap.begin(), heap.end(), Cmp());
  while (!heap.empty()) {
    std::pop_heap(heap.begin(), heap.end(), Cmp());
    Node n = heap.back(); heap.pop_back();
    out.push_back(n.v);
    int r = n.r;
    runs[r].idx++;
    if (runs[r].idx < runs[r].end) {
      heap.push_back({in[runs[r].idx], r});
      std::push_heap(heap.begin(), heap.end(), Cmp());
    }
  }
}

// ----------------------------- Sample Sort -----------------------------

template<class T>
static void sample_sort(std::vector<T>& local, MPI_Comm comm) {
  int rank, P; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&P);
  const bool is_double = std::is_same<T,double>::value;

  // ---- Local sort (comp_large) ----
  CALI_MARK_BEGIN("comp");
  CALI_MARK_BEGIN("comp_large");
  std::sort(local.begin(), local.end());
  CALI_MARK_END("comp_large");
  CALI_MARK_END("comp");

  // ---- Local sampling (comp_small) ----
  CALI_MARK_BEGIN("comp");
  CALI_MARK_BEGIN("comp_small");
  int s = std::max(1, P);
  std::vector<T> samples;
  samples.reserve((size_t)std::min<long long>(local.size(), (long long)s));
  if (!local.empty()) {
    for (int i=0;i<s;++i) {
      size_t idx = (size_t)((i+1) * (local.size() / (double)(s+1)));
      if (idx >= local.size()) idx = local.size()-1;
      samples.push_back(local[idx]);
    }
  }
  CALI_MARK_END("comp_small");
  CALI_MARK_END("comp");

  // ---- Gather samples (comm_small) ----
  CALI_MARK_BEGIN("comm");
  CALI_MARK_BEGIN("comm_small");
  int loc_s = (int)samples.size();
  std::vector<int> all_counts(P);
  MPI_Allgather(&loc_s, 1, MPI_INT, all_counts.data(), 1, MPI_INT, comm);
  int total_s = std::accumulate(all_counts.begin(), all_counts.end(), 0);
  std::vector<int> displs(P,0);
  for (int i=1;i<P;++i) displs[i] = displs[i-1] + all_counts[i-1];

  std::vector<T> all_samples((size_t)total_s);
  MPI_Datatype Ttype = is_double ? MPI_DOUBLE : MPI_LONG_LONG;
  MPI_Allgatherv(samples.data(), loc_s, Ttype,
                 all_samples.data(), all_counts.data(), displs.data(), Ttype,
                 comm);
  CALI_MARK_END("comm_small");
  CALI_MARK_END("comm");

  // ---- Sort samples and choose splitters (comp_small) ----
  CALI_MARK_BEGIN("comp");
  CALI_MARK_BEGIN("comp_small");
  std::sort(all_samples.begin(), all_samples.end());
  std::vector<T> splitters;
  splitters.reserve((size_t)std::max(0, P-1));
  if (!all_samples.empty()) {
    for (int i=1;i<P;++i) {
      size_t idx = (size_t) std::llround((long double)i * (long double)all_samples.size() / (long double)P);
      if (idx>=all_samples.size()) idx = all_samples.size()-1;
      splitters.push_back(all_samples[idx]);
    }
  } else {
    splitters.assign((size_t)std::max(0,P-1), T());
  }
  CALI_MARK_END("comp_small");
  CALI_MARK_END("comp");

  // ---- Partition local array into P buckets using splitters (comp_small) ----
  CALI_MARK_BEGIN("comp");
  CALI_MARK_BEGIN("comp_small");
  std::vector<int> send_counts(P,0);
  if (!local.empty()) {
    std::vector<size_t> bounds; bounds.reserve((size_t)P+1);
    bounds.push_back(0);
    for (int i=0;i<P-1;++i) {
      typename std::vector<T>::iterator it = std::upper_bound(local.begin()+bounds.back(), local.end(), splitters[i]);
      bounds.push_back((size_t)std::distance(local.begin(), it));
    }
    bounds.push_back(local.size());
    for (int i=0;i<P;++i) {
      send_counts[i] = (int)(bounds[i+1]-bounds[i]);
    }
  } else {
    std::fill(send_counts.begin(), send_counts.end(), 0);
  }

  std::vector<int> send_displs(P,0);
  for (int i=1;i<P;++i) send_displs[i] = send_displs[i-1] + send_counts[i-1];
  std::vector<T> sendbuf(local.size());
  {
    if (!local.empty()) {
      std::vector<size_t> bounds; bounds.reserve((size_t)P+1);
      bounds.push_back(0);
      for (int i=0;i<P-1;++i) {
        typename std::vector<T>::iterator it = std::upper_bound(local.begin()+bounds.back(), local.end(), splitters[i]);
        bounds.push_back((size_t)std::distance(local.begin(), it));
      }
      bounds.push_back(local.size());
      for (int i=0;i<P;++i) {
        size_t b = bounds[i], e = bounds[i+1];
        std::copy(local.begin()+b, local.begin()+e, sendbuf.begin()+send_displs[i]);
      }
    }
  }
  CALI_MARK_END("comp_small");
  CALI_MARK_END("comp");

  // ---- All-to-all sizes (comm_small) ----
  CALI_MARK_BEGIN("comm");
  CALI_MARK_BEGIN("comm_small");
  std::vector<int> recv_counts(P,0);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT,
               recv_counts.data(), 1, MPI_INT, comm);
  std::vector<int> recv_displs(P,0);
  int recv_total = 0;
  for (int i=0;i<P;++i){ recv_displs[i]=recv_total; recv_total+=recv_counts[i]; }
  std::vector<T> recvbuf((size_t)recv_total);
  CALI_MARK_END("comm_small");
  CALI_MARK_END("comm");

  // ---- All-to-allv exchange (comm_large) ----
  CALI_MARK_BEGIN("comm");
  CALI_MARK_BEGIN("comm_large");
  MPI_Alltoallv(sendbuf.data(), send_counts.data(), send_displs.data(), Ttype,
                recvbuf.data(), recv_counts.data(), recv_displs.data(), Ttype,
                comm);
  CALI_MARK_END("comm_large");
  CALI_MARK_END("comm");

  // ---- Final local k-way merge (comp_large) ----
  CALI_MARK_BEGIN("comp");
  CALI_MARK_BEGIN("comp_large");
  std::vector<T> final_sorted; final_sorted.reserve(recvbuf.size());
  kway_merge(recvbuf, recv_counts, final_sorted);
  local.swap(final_sorted);
  CALI_MARK_END("comp_large");
  CALI_MARK_END("comp");
}

// ----------------------------- Data init helpers (tag dispatch) -----------------------------

template<typename T> struct DataInit;

template<> struct DataInit<int64_t> {
  static void fill(std::vector<int64_t>& v, const Options& opt, int rank) {
    std::mt19937_64 gen(opt.seed + (uint64_t)rank);
    std::uniform_int_distribution<int64_t> dist(-1000000000LL, 1000000000LL);
    if (opt.input=="random") {
      for (size_t i=0;i<v.size();++i) v[i] = dist(gen);
    } else if (opt.input=="sorted") {
      for (size_t i=0;i<v.size();++i) v[i] = (int64_t)(rank) * 1000000000LL + (int64_t)i;
    } else if (opt.input=="reversed") {
      for (size_t i=0;i<v.size();++i) v[i] = -(int64_t)((rank) * 1000000000LL + (int64_t)i);
    } else if (opt.input=="perturbed1") {
      for (size_t i=0;i<v.size();++i) v[i] = (int64_t)i;
      size_t swaps = v.size()/100 + 1;
      for (size_t s=0;s<swaps;++s) { size_t a = gen()%v.size(), b = gen()%v.size(); std::swap(v[a], v[b]); }
    } else {
      for (size_t i=0;i<v.size();++i) v[i] = dist(gen);
    }
  }
};

template<> struct DataInit<double> {
  static void fill(std::vector<double>& v, const Options& opt, int rank) {
    std::mt19937_64 gen(opt.seed + (uint64_t)rank);
    std::uniform_real_distribution<double> dist(-1e6, 1e6);
    if (opt.input=="random") {
      for (size_t i=0;i<v.size();++i) v[i] = dist(gen);
    } else if (opt.input=="sorted") {
      for (size_t i=0;i<v.size();++i) v[i] = (double)rank * 1e6 + (double)i;
    } else if (opt.input=="reversed") {
      for (size_t i=0;i<v.size();++i) v[i] = -(double)((double)rank * 1e6 + (double)i);
    } else if (opt.input=="perturbed1") {
      for (size_t i=0;i<v.size();++i) v[i] = (double)i;
      size_t swaps = v.size()/100 + 1;
      for (size_t s=0;s<swaps;++s) { size_t a = gen()%v.size(), b = gen()%v.size(); std::swap(v[a], v[b]); }
    } else {
      for (size_t i=0;i<v.size();++i) v[i] = dist(gen);
    }
  }
};

template<class T>
static std::vector<T> make_local_data(const Options& opt, MPI_Comm comm) {
  int rank, P; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&P);
  long long Nloc_base = opt.N / P;
  int remainder = (int)(opt.N % P);
  long long nloc = Nloc_base + (rank < remainder ? 1 : 0);

  std::vector<T> v; v.resize((size_t)nloc);
  CALI_MARK_BEGIN("data_init_runtime");
  DataInit<T>::fill(v, opt, rank);
  CALI_MARK_END("data_init_runtime");
  return v;
}

// ----------------------------- Main -----------------------------

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank, P; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&P);

  CALI_CXX_MARK_FUNCTION; // "main"

  cali::ConfigManager mgr;
  mgr.start();

  // Adiak metadata
  adiak::init(nullptr);
  Options opt = parse_cli(argc, argv, rank);
  adiak::launchdate();
  adiak::libraries();
  adiak::cmdline();
  adiak::clustername();
  adiak::value("algorithm", "sample-sort");
  adiak::value("programming_model", "mpi");
  adiak::value("data_type", opt.datatype);
  adiak::value("size_of_data_type", (int)(opt.datatype=="double" ? sizeof(double) : sizeof(int64_t)));
  adiak::value("input_size", (long long)opt.N);
  adiak::value("input_type", opt.input);
  adiak::value("num_procs", P);
  adiak::value("scalability", opt.scalability);
  adiak::value("group_num", opt.group_num);
  adiak::value("implementation_source", "handwritten");

  // Choose datatype
  bool use_double = (opt.datatype=="double");

  // World communicator
  MPI_Comm comm = MPI_COMM_WORLD;

  // Generate data
  std::vector<int64_t> local_i64;
  std::vector<double>  local_d;
  if (!use_double) local_i64 = make_local_data<int64_t>(opt, comm);
  else             local_d   = make_local_data<double>(opt, comm);

  // Run sort
  if (!use_double) sample_sort<int64_t>(local_i64, comm);
  else             sample_sort<double>(local_d, comm);

  // Correctness check (global)
  CALI_MARK_BEGIN("correctness_check");
  bool ok = false;
  if (!use_double) ok = global_sorted_check<int64_t>(local_i64, comm);
  else             ok = global_sorted_check<double>(local_d, comm);
  CALI_MARK_END("correctness_check");

  if (!rank) {
    std::cout << "Global sorted check: " << (ok ? "PASSED" : "FAILED") << std::endl;
  }

  mgr.flush();
  MPI_Finalize();
  return ok ? 0 : 2;
}
