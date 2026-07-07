#ifndef PVFMM_WRAP_H
#define PVFMM_WRAP_H

// ---------------------------------------------------------------------------
// Single-file PVFMM wrapper.
//
// This header carries BOTH the C-callable declarations (extern "C") AND their
// C++ definitions. It is injected into the mwrap-generated gateway (pvfmm.c)
// via the makefile's `-include pvfmm_wrap.h`, so the gateway becomes the one
// and only translation unit that defines these symbols -> no separate
// pvfmm_wrap.cpp, no duplicate symbols.
//
// REQUIREMENT: the gateway must be compiled as C++ (e.g. g++ / g++-14), because
// the bodies use std::vector, PVFMM templates, MPI, etc. The makefile compiles
// pvfmm.c as C++ on BOTH the MATLAB and the Octave paths. If you ever revert a
// build rule to plain C (e.g. `gcc -std=c99`) the #error below fires.
// ---------------------------------------------------------------------------

#ifndef __cplusplus
#error "pvfmm_wrap.h contains the C++ PVFMM implementation; compile the mwrap gateway (pvfmm.c) as C++ (e.g. with g++), not as C (e.g. gcc -std=c99)."
#endif

#include <stdint.h>
#include <mpi.h>
#include <omp.h>
#include <pvfmm.hpp>
#include <sctl.hpp>          // SCTL ParticleFMM + Stokes3D_* kernels (pvfmm.hpp does not pull these)
#include <vector>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using vec = std::vector<double>;

struct MPI_Context {
  MPI_Comm comm;
  int rank;
  int size;

  MPI_Context(MPI_Comm c, int r, int n) : comm(c), rank(r), size(n) {}
};

// Internal C++ helpers (defined below; called only within this TU).
void laplace_slp_naive(vec& sl_coord, vec& sl_den, vec& trg_coord, vec& trg_value, MPI_Comm& comm);
void laplace_slp_fmm( const vec&  sl_coord, const vec& sl_den, const vec& trg_coord, vec& trg_value, const int mult_order, const size_t max_pts, const size_t mem_size, MPI_Comm& comm);
void laplace_slpn_fmm( const vec&  sl_coord, const vec& sl_den, const vec& trg_coord, const vec& trg_normal, vec& trg_value, const int mult_order, const size_t max_pts, const size_t mem_size, MPI_Comm& comm);
void stokes_slp_fmm( const vec&  sl_coord, const vec& sl_den, const vec& trg_coord, vec& trg_value, const int mult_order, const size_t max_pts, const size_t mem_size, MPI_Comm& comm);
void stokes_slpn_fmm( const vec&  sl_coord, const vec& sl_den, const vec& trg_coord, const vec& trg_normal, vec& trg_value, const int mult_order, const size_t max_pts, const size_t mem_size, MPI_Comm& comm);
void stokes_dlp_fmm( const vec&  dl_coord, const vec& dl_sigma, const vec& dl_normal, const vec& trg_coord, vec& trg_value, const int mult_order, const size_t max_pts, const size_t mem_size, MPI_Comm& comm);
void stokes_dlpn_fmm( const vec&  dl_coord, const vec& dl_sigma, const vec& dl_normal, const vec& trg_coord, const vec& trg_normal, vec& trg_value, const int mult_order, const size_t max_pts, const size_t mem_size, MPI_Comm& comm);
void laplace_dlp_fmm( const vec&  dl_coord, const vec& dl_normal, const vec& dl_sigma, const vec& trg_coord, vec& trg_value, const int mult_order, const size_t max_pts, const size_t mem_size, MPI_Comm& comm);
void laplace_dlpn_fmm( const vec&  dl_coord, const vec& dl_normal, const vec& dl_sigma, const vec& trg_coord, const vec& trg_normal, vec& trg_value, const int mult_order, const size_t max_pts, const size_t mem_size, MPI_Comm& comm);
void helmholtz_slp_fmm( const vec&  sl_coord, const vec& sl_den, const vec& trg_coord, vec& trg_value, const int mult_order, const size_t max_pts, const size_t mem_size, MPI_Comm& comm);
extern "C" void mpi_init_wrap(const char* args_string, uint64_t *handle, int* flag) {
  int initialized = 0;
  MPI_Initialized(&initialized);

  int rank, size;

  // --- 1. args_string to argc/argv ---
  std::vector<std::string> tokens;
  std::istringstream iss(args_string ? args_string : "MATLAB_MPI");
  std::string t;
  while (iss >> t) tokens.push_back(t);

  for (size_t i = 0; i < tokens.size(); ++i) {
    if (tokens[i] == "-omp" && (i + 1) < tokens.size()) {
      int num_threads = std::atoi(tokens[i + 1].c_str());
      if (num_threads > 0) {
        omp_set_num_threads(num_threads);
        std::printf("[MATLAB_MPI INFO] OpenMP threads set to %d\n", num_threads);
      }
    }
  }

  if (!initialized) {

    int argc = static_cast<int>(tokens.size());
    std::vector<char*> argv_vec;
    for (auto& s : tokens) {
      argv_vec.push_back(strdup(s.c_str()));
    }
    argv_vec.push_back(nullptr);
    char** argv = argv_vec.data();

    // --- 2. MPI_Init ---
    *flag = MPI_Init(&argc, &argv);

    // strdup, not sure
    for (auto p : argv_vec) if(p) free(p);

    // --- 3. create handle ---
    if (*flag == MPI_SUCCESS) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      MPI_Context* ctx = new MPI_Context(MPI_COMM_WORLD, rank, size);
      *handle = reinterpret_cast<uint64_t>(ctx);
      std::printf("[MATLAB_MPI INFO] Initialized. Rank %d/%d\n", ctx->rank, ctx->size);
    }
  } else {
    if (handle && *handle != 0) {
      *flag = MPI_SUCCESS;
      std::printf("[MATLAB_MPI INFO] Handle already exists. Skipping allocation.\n");
      return;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Context* ctx = new MPI_Context(MPI_COMM_WORLD, rank, size);
    *handle = reinterpret_cast<uint64_t>(ctx);
    *flag = MPI_SUCCESS;
    std::printf("[MATLAB_MPI INFO] MPI already active. Restored Context for Rank %d/%d\n", ctx->rank, ctx->size);
    // if already created... handle
  }
}
// // Simple version...
// extern "C" void mpi_init_wrap(int* flag) {
//   int initialized = 0;
//   MPI_Initialized(&initialized);
//   if (!initialized) {
//     *flag = MPI_Init(NULL, NULL);
//     return;
//   }
//   *flag = MPI_SUCCESS;
// }

extern "C" void mpi_finalize_wrap(uint64_t *handle, int* flag) {
  // 0. inout flag may come from user request
  int user_request = (flag != nullptr) ? *flag : 0;

  // 1. clear memory (Handle)
  if (handle && *handle != 0) {
    MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
    delete ctx;
    *handle = 0;
  }
  // 2. clear MPI
  int initialized = 0;
  int finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (!initialized || finalized) {
    *flag = MPI_SUCCESS;
    // print notice for Matlab users
    std::printf("\n****************************************************************\n");
    std::printf("[MATLAB_MPI INFO] MPI is not active or already finalized.\n");
    std::printf("****************************************************************\n\n");
    fflush(stdout);
    return;
  }

  if (user_request == 99) {
    std::printf("\n****************************************************************\n");
    std::printf("[MATLAB_MPI INFO] CRITICAL SHUTDOWN (Flag 99) detected.\n");
    std::printf("Executing MPI_Finalize(). Expecting MATLAB thread instability...\n");
    std::printf("To run MPI tasks again or MATLAB crashes now, please RESTART another MATLAB session.\n");
    std::printf("****************************************************************\n\n");
    fflush(stdout);
    *flag = MPI_Finalize();
    std::exit(0);
  } else {
    *flag = MPI_SUCCESS;
    std::printf("\n****************************************************************\n");
    std::printf("[MATLAB_MPI INFO] Context memory has been cleared.\n");
    std::printf("[KEEPING MPI ACTIVE] Skipping MPI_Finalize() for stability.\n");
    std::printf("To TOTALLY shutdown MPI (and risk a crash), call with flag = 99:\n");
    std::printf("Note: After a complete MPI shutdown, a full restart of MATLAB is required to re-initialize the MPI environment.\n");
    std::printf("   >> handle = 0; flag = 99; [handle, flag] = MPI_Finalize_mex(handle, flag);\n");
    std::printf("****************************************************************\n\n");
    fflush(stdout);
  }

}
// // Simple version...
// extern "C" void mpi_finalize_wrap(int* flag) {
//   int finalized = 0;
//   MPI_Finalized(&finalized);
//   if (!finalized) {
//     *flag = MPI_Finalize();
//     return;
//   }
//   *flag = MPI_SUCCESS;
// }

extern "C" void mpi_comm_rank_wrap(int* rank) {
  if (!rank) {
    return;
  }
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, rank);
}

extern "C" void mpi_comm_size_wrap(int* size) {
  if (!size) {
    return;
  }
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, size);
}

extern "C" void laplace_slp_naive_wrap(uint64_t *handle,
                                    const double* sl_coord, const double* sl_den,
                                    const double* trg_coord, int* n_sl, int* n_trg,
                                    double* trg_value) {
  //
  if (!handle || *handle == 0 ||!sl_coord || !sl_den || !trg_coord || !trg_value || !n_sl || !n_trg) {
    return;
  }
  if (*n_sl <= 0 || *n_trg <= 0) {
    return;
  }
  //
  MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
  MPI_Comm comm = ctx->comm;
  //
  const size_t sl_count = static_cast<size_t>(*n_sl) * PVFMM_COORD_DIM;
  const size_t trg_count = static_cast<size_t>(*n_trg) * PVFMM_COORD_DIM;
  vec sl_coord_v(sl_coord, sl_coord + sl_count);
  vec sl_den_v(sl_den, sl_den + static_cast<size_t>(*n_sl));
  vec trg_coord_v(trg_coord, trg_coord + trg_count);
  vec trg_value_v;
  //
  laplace_slp_naive(sl_coord_v, sl_den_v, trg_coord_v, trg_value_v, comm);
  if (trg_value_v.size() >= static_cast<size_t>(*n_trg)) {
    std::copy(trg_value_v.begin(), trg_value_v.begin() + *n_trg, trg_value);
  }
}

extern "C" void laplace_slp_fmm_wrap(uint64_t *handle,
                                    const double* sl_coord, const double* sl_den,
                                    const double* trg_coord, int* n_sl, int* n_trg,
                                    int* mult_order, uint64_t *mem_size_in, int* max_pts_in,
                                    double* trg_value) {
  //
  if (!handle || *handle == 0 ||!sl_coord || !sl_den || !trg_coord || !trg_value || !n_sl || !n_trg) {
    return;
  }
  if (*n_sl <= 0 || *n_trg <= 0) {
    return;
  }
  //
  MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
  MPI_Comm comm = ctx->comm;
  //
  size_t mem_size = static_cast<size_t>(*mem_size_in);
  size_t max_pts  = static_cast<size_t>(*max_pts_in);
  //
  const size_t sl_den_count = static_cast<size_t>(*n_sl);
  const size_t trg_value_count = static_cast<size_t>(*n_trg);
  const size_t sl_count     = sl_den_count * PVFMM_COORD_DIM;
  const size_t trg_count    = trg_value_count * PVFMM_COORD_DIM;
  vec sl_coord_v(sl_coord, sl_coord + sl_count);
  vec sl_den_v(sl_den, sl_den + sl_den_count);
  vec trg_coord_v(trg_coord, trg_coord + trg_count);
  vec trg_value_v;
  //
  laplace_slp_fmm(sl_coord_v, sl_den_v, trg_coord_v, trg_value_v,
                  *mult_order, max_pts, mem_size, comm);
  if (trg_value_v.size() >= trg_value_count) {
    std::copy(trg_value_v.begin(), trg_value_v.begin() + trg_value_count, trg_value);
  }
}

extern "C" void stokes_slp_fmm_wrap(uint64_t *handle,
                                    const double* sl_coord, const double* sl_den,
                                    const double* trg_coord, int* n_sl, int* n_trg,
                                    int* mult_order, uint64_t *mem_size_in, int* max_pts_in,
                                    double* trg_value) {
  //
  if (!handle || *handle == 0 ||!sl_coord || !sl_den || !trg_coord || !trg_value || !n_sl || !n_trg) {
    return;
  }
  if (*n_sl <= 0 || *n_trg <= 0) {
    return;
  }
  //
  MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
  MPI_Comm comm = ctx->comm;
  //
  size_t mem_size = static_cast<size_t>(*mem_size_in);
  size_t max_pts  = static_cast<size_t>(*max_pts_in);
  //
  const size_t sl_den_count = static_cast<size_t>(*n_sl) * PVFMM_COORD_DIM;
  const size_t trg_value_count = static_cast<size_t>(*n_trg) * PVFMM_COORD_DIM;
  const size_t sl_count     = static_cast<size_t>(*n_sl) * PVFMM_COORD_DIM;
  const size_t trg_count    = static_cast<size_t>(*n_trg) * PVFMM_COORD_DIM;
  vec sl_coord_v(sl_coord, sl_coord + sl_count);
  vec sl_den_v(sl_den, sl_den + sl_den_count);
  vec trg_coord_v(trg_coord, trg_coord + trg_count);
  vec trg_value_v;
  //
  stokes_slp_fmm(sl_coord_v, sl_den_v, trg_coord_v, trg_value_v,
                 *mult_order, max_pts, mem_size, comm);
  if (trg_value_v.size() >= trg_value_count) {
    std::copy(trg_value_v.begin(), trg_value_v.begin() + trg_value_count, trg_value);
  }
}

extern "C" void stokes_slpn_fmm_wrap(uint64_t *handle,
                                     const double* sl_coord, const double* sl_den,
                                     const double* trg_coord, const double* trg_normal,
                                     int* n_sl, int* n_trg,
                                     int* mult_order, uint64_t *mem_size_in, int* max_pts_in,
                                     double* trg_value) {
  if (!handle || *handle == 0 || !sl_coord || !sl_den || !trg_coord || !trg_normal ||
      !trg_value || !n_sl || !n_trg) {
    return;
  }
  if (*n_sl <= 0 || *n_trg <= 0) {
    return;
  }
  MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
  MPI_Comm comm = ctx->comm;
  size_t mem_size = static_cast<size_t>(*mem_size_in);
  size_t max_pts  = static_cast<size_t>(*max_pts_in);
  const size_t sl_count        = static_cast<size_t>(*n_sl)  * PVFMM_COORD_DIM;   // 3-vec force
  const size_t trg_count       = static_cast<size_t>(*n_trg) * PVFMM_COORD_DIM;   // coords / normals
  const size_t trg_value_count = static_cast<size_t>(*n_trg) * PVFMM_COORD_DIM;   // 3-vec traction
  vec sl_coord_v(sl_coord, sl_coord + sl_count);
  vec sl_den_v(sl_den, sl_den + sl_count);
  vec trg_coord_v(trg_coord, trg_coord + trg_count);
  vec trg_normal_v(trg_normal, trg_normal + trg_count);
  vec trg_value_v;
  stokes_slpn_fmm(sl_coord_v, sl_den_v, trg_coord_v, trg_normal_v, trg_value_v,
                  *mult_order, max_pts, mem_size, comm);
  if (trg_value_v.size() >= trg_value_count) {
    std::copy(trg_value_v.begin(), trg_value_v.begin() + trg_value_count, trg_value);
  }
}

extern "C" void stokes_dlp_fmm_wrap(uint64_t *handle,
                                    const double* dl_coord, const double* dl_sigma,
                                    const double* dl_normal, const double* trg_coord,
                                    int* n_dl, int* n_trg,
                                    int* mult_order, uint64_t *mem_size_in, int* max_pts_in,
                                    double* trg_value) {
  if (!handle || *handle == 0 || !dl_coord || !dl_sigma || !dl_normal || !trg_coord ||
      !trg_value || !n_dl || !n_trg) {
    return;
  }
  if (*n_dl <= 0 || *n_trg <= 0) {
    return;
  }
  MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
  MPI_Comm comm = ctx->comm;
  size_t mem_size = static_cast<size_t>(*mem_size_in);
  size_t max_pts  = static_cast<size_t>(*max_pts_in);
  const size_t dl_count        = static_cast<size_t>(*n_dl)  * PVFMM_COORD_DIM;   // coords / sigma / normals
  const size_t trg_count       = static_cast<size_t>(*n_trg) * PVFMM_COORD_DIM;   // coords
  const size_t trg_value_count = static_cast<size_t>(*n_trg) * PVFMM_COORD_DIM;   // 3-vec velocity
  vec dl_coord_v(dl_coord, dl_coord + dl_count);
  vec dl_sigma_v(dl_sigma, dl_sigma + dl_count);
  vec dl_normal_v(dl_normal, dl_normal + dl_count);
  vec trg_coord_v(trg_coord, trg_coord + trg_count);
  vec trg_value_v;
  stokes_dlp_fmm(dl_coord_v, dl_sigma_v, dl_normal_v, trg_coord_v, trg_value_v,
                 *mult_order, max_pts, mem_size, comm);
  if (trg_value_v.size() >= trg_value_count) {
    std::copy(trg_value_v.begin(), trg_value_v.begin() + trg_value_count, trg_value);
  }
}

extern "C" void stokes_dlpn_fmm_wrap(uint64_t *handle,
                                     const double* dl_coord, const double* dl_sigma,
                                     const double* dl_normal, const double* trg_coord,
                                     const double* trg_normal, int* n_dl, int* n_trg,
                                     int* mult_order, uint64_t *mem_size_in, int* max_pts_in,
                                     double* trg_value) {
  if (!handle || *handle == 0 || !dl_coord || !dl_sigma || !dl_normal || !trg_coord ||
      !trg_normal || !trg_value || !n_dl || !n_trg) {
    return;
  }
  if (*n_dl <= 0 || *n_trg <= 0) {
    return;
  }
  MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
  MPI_Comm comm = ctx->comm;
  size_t mem_size = static_cast<size_t>(*mem_size_in);
  size_t max_pts  = static_cast<size_t>(*max_pts_in);
  const size_t dl_count        = static_cast<size_t>(*n_dl)  * PVFMM_COORD_DIM;   // coords / sigma / normals
  const size_t trg_count       = static_cast<size_t>(*n_trg) * PVFMM_COORD_DIM;   // coords / normals
  const size_t trg_value_count = static_cast<size_t>(*n_trg) * PVFMM_COORD_DIM;   // 3-vec traction
  vec dl_coord_v(dl_coord, dl_coord + dl_count);
  vec dl_sigma_v(dl_sigma, dl_sigma + dl_count);
  vec dl_normal_v(dl_normal, dl_normal + dl_count);
  vec trg_coord_v(trg_coord, trg_coord + trg_count);
  vec trg_normal_v(trg_normal, trg_normal + trg_count);
  vec trg_value_v;
  stokes_dlpn_fmm(dl_coord_v, dl_sigma_v, dl_normal_v, trg_coord_v, trg_normal_v, trg_value_v,
                  *mult_order, max_pts, mem_size, comm);
  if (trg_value_v.size() >= trg_value_count) {
    std::copy(trg_value_v.begin(), trg_value_v.begin() + trg_value_count, trg_value);
  }
}

extern "C" void laplace_slpn_fmm_wrap(uint64_t *handle,
                                    const double* sl_coord, const double* sl_den,
                                    const double* trg_coord, const double* trg_normal,
                                    int* n_sl, int* n_trg,
                                    int* mult_order, uint64_t *mem_size_in, int* max_pts_in,
                                    double* trg_value) {
  if (!handle || *handle == 0 || !sl_coord || !sl_den || !trg_coord || !trg_normal ||
      !trg_value || !n_sl || !n_trg) {
    return;
  }
  if (*n_sl <= 0 || *n_trg <= 0) {
    return;
  }
  //
  MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
  MPI_Comm comm = ctx->comm;
  //
  size_t mem_size = static_cast<size_t>(*mem_size_in);
  size_t max_pts  = static_cast<size_t>(*max_pts_in);
  //
  const size_t sl_den_count = static_cast<size_t>(*n_sl);
  const size_t trg_value_count = static_cast<size_t>(*n_trg);
  const size_t sl_count     = sl_den_count * PVFMM_COORD_DIM;
  const size_t trg_count    = trg_value_count * PVFMM_COORD_DIM;
  vec sl_coord_v(sl_coord, sl_coord + sl_count);
  vec sl_den_v(sl_den, sl_den + sl_den_count);
  vec trg_coord_v(trg_coord, trg_coord + trg_count);
  vec trg_normal_v(trg_normal, trg_normal + trg_count);
  vec trg_value_v;
  //
  laplace_slpn_fmm(sl_coord_v, sl_den_v, trg_coord_v, trg_normal_v, trg_value_v,
                   *mult_order, max_pts, mem_size, comm);
  if (trg_value_v.size() >= trg_value_count) {
    std::copy(trg_value_v.begin(), trg_value_v.begin() + trg_value_count, trg_value);
  }
}

extern "C" void laplace_dlp_fmm_wrap(uint64_t *handle,
                                    const double* dl_coord, const double* dl_normal, const double* dl_sigma,
                                    const double* trg_coord, int* n_dl, int* n_trg,
                                    int* mult_order, uint64_t *mem_size_in, int* max_pts_in,
                                    double* trg_value) {
  //
  if (!handle || *handle == 0 ||!dl_coord ||!dl_normal || !dl_sigma || !trg_coord || !trg_value || !n_dl || !n_trg) {
    return;
  }
  if (*n_dl <= 0 || *n_trg <= 0) {
    return;
  }
  //
  MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
  MPI_Comm comm = ctx->comm;
  //
  size_t mem_size = static_cast<size_t>(*mem_size_in);
  size_t max_pts  = static_cast<size_t>(*max_pts_in);
  //
  const size_t dl_sigma_count = static_cast<size_t>(*n_dl);
  const size_t trg_value_count = static_cast<size_t>(*n_trg);
  const size_t dl_count = dl_sigma_count * PVFMM_COORD_DIM;
  const size_t trg_count = trg_value_count * PVFMM_COORD_DIM;
  vec dl_coord_v(dl_coord, dl_coord + dl_count);
  vec dl_normal_v(dl_normal, dl_normal + dl_count);
  vec dl_sigma_v(dl_sigma, dl_sigma + dl_sigma_count);
  vec trg_coord_v(trg_coord, trg_coord + trg_count);
  vec trg_value_v;
  //
  laplace_dlp_fmm(dl_coord_v, dl_normal_v, dl_sigma_v, trg_coord_v, trg_value_v,
                  *mult_order, max_pts, mem_size, comm);
  if (trg_value_v.size() >= trg_value_count) {
    std::copy(trg_value_v.begin(), trg_value_v.begin() + trg_value_count, trg_value);
  }
}

extern "C" void laplace_dlpn_fmm_wrap(uint64_t *handle,
                                    const double* dl_coord, const double* dl_normal, const double* dl_sigma,
                                    const double* trg_coord, const double* trg_normal, int* n_dl, int* n_trg,
                                    int* mult_order, uint64_t *mem_size_in, int* max_pts_in,
                                    double* trg_value) {
  if (!handle || *handle == 0 || !dl_coord || !dl_normal || !dl_sigma || !trg_coord ||
      !trg_normal || !trg_value || !n_dl || !n_trg) {
    return;
  }
  if (*n_dl <= 0 || *n_trg <= 0) {
    return;
  }
  MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
  MPI_Comm comm = ctx->comm;
  size_t mem_size = static_cast<size_t>(*mem_size_in);
  size_t max_pts  = static_cast<size_t>(*max_pts_in);
  const size_t dl_sigma_count  = static_cast<size_t>(*n_dl);
  const size_t trg_value_count = static_cast<size_t>(*n_trg);          // scalar DLPn per target
  const size_t dl_count  = dl_sigma_count  * PVFMM_COORD_DIM;
  const size_t trg_count = trg_value_count * PVFMM_COORD_DIM;
  vec dl_coord_v(dl_coord, dl_coord + dl_count);
  vec dl_normal_v(dl_normal, dl_normal + dl_count);
  vec dl_sigma_v(dl_sigma, dl_sigma + dl_sigma_count);
  vec trg_coord_v(trg_coord, trg_coord + trg_count);
  vec trg_normal_v(trg_normal, trg_normal + trg_count);
  vec trg_value_v;
  laplace_dlpn_fmm(dl_coord_v, dl_normal_v, dl_sigma_v, trg_coord_v, trg_normal_v, trg_value_v,
                   *mult_order, max_pts, mem_size, comm);
  if (trg_value_v.size() >= trg_value_count) {
    std::copy(trg_value_v.begin(), trg_value_v.begin() + trg_value_count, trg_value);
  }
}

extern "C" void helmholtz_slp_fmm_wrap(uint64_t *handle,
                                    const double* sl_coord, const double* sl_den,
                                    const double* trg_coord, int* n_sl, int* n_trg,
                                    int* mult_order, uint64_t *mem_size_in, int* max_pts_in,
                                    double* trg_value) {
  if (!handle || *handle == 0 || !sl_coord || !sl_den || !trg_coord || !trg_value ||
      !n_sl || !n_trg) {
    return;
  }
  if (*n_sl <= 0 || *n_trg <= 0) {
    return;
  }
  //
  MPI_Context* ctx = reinterpret_cast<MPI_Context*>(*handle);
  MPI_Comm comm = ctx->comm;
  //
  size_t mem_size = static_cast<size_t>(*mem_size_in);
  size_t max_pts  = static_cast<size_t>(*max_pts_in);
  //
  const size_t sl_den_count = static_cast<size_t>(*n_sl) * 2;
  const size_t trg_value_count = static_cast<size_t>(*n_trg) * 2;
  const size_t sl_count     = static_cast<size_t>(*n_sl) * PVFMM_COORD_DIM;
  const size_t trg_count    = static_cast<size_t>(*n_trg) * PVFMM_COORD_DIM;
  vec sl_coord_v(sl_coord, sl_coord + sl_count);
  vec sl_den_v(sl_den, sl_den + sl_den_count);
  vec trg_coord_v(trg_coord, trg_coord + trg_count);
  vec trg_value_v;
  //
  helmholtz_slp_fmm(sl_coord_v, sl_den_v, trg_coord_v, trg_value_v,
                    *mult_order, max_pts, mem_size, comm);
  if (trg_value_v.size() >= trg_value_count) {
    std::copy(trg_value_v.begin(), trg_value_v.begin() + trg_value_count, trg_value);
  }
}

//
void laplace_slp_naive(vec& sl_coord, vec& sl_den, vec& trg_coord, vec& trg_value,
                          MPI_Comm& comm){
  (void)comm;
  const size_t Ns = sl_coord.size() / PVFMM_COORD_DIM;
  const size_t Nt = trg_coord.size() / PVFMM_COORD_DIM;
  trg_value.assign(Nt, 0.0);

  const double oofp = 1.0 / (4.0 * M_PI);
  const double tol2 = 1e-28;

  // #pragma omp parallel for schedule(static)
  for (size_t t = 0; t < Nt; ++t) {
    double acc = 0.0;
    const double xt0 = trg_coord[t * 3 + 0];
    const double xt1 = trg_coord[t * 3 + 1];
    const double xt2 = trg_coord[t * 3 + 2];

    for (size_t s = 0; s < Ns; ++s) {
      const double dx0 = xt0 - sl_coord[s * 3 + 0];
      const double dx1 = xt1 - sl_coord[s * 3 + 1];
      const double dx2 = xt2 - sl_coord[s * 3 + 2];
      const double r2 = dx0 * dx0 + dx1 * dx1 + dx2 * dx2;
      if (r2 > tol2) {
        acc += sl_den[s] / std::sqrt(r2);
      }
    }
    trg_value[t] = acc * oofp;
  }
}


// kernel_fn.ker_poten --> SLP in kernel.txx
void laplace_slp_fmm( const vec&  sl_coord, const vec& sl_den,
                      const vec& trg_coord, vec& trg_value,
                      const int mult_order, const size_t max_pts, const size_t mem_size,
                      MPI_Comm& comm){
  const pvfmm::Kernel<double>& kernel_fn = pvfmm::LaplaceKernel<double>::potential();
  const size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  //
  sctl::Comm sctl_comm(comm);
  // size_t max_pts = 600;
  // Empty double-layer inputs for SLP-only run.
  vec dl_coord, dl_den;
  auto* tree = PtFMM_CreateTree(sl_coord, sl_den, dl_coord, dl_den,
                                trg_coord, sctl_comm, (int)max_pts, pvfmm::FreeSpace);
  //
  pvfmm::PtFMM<double> matrices;
  matrices.Initialize(mult_order, sctl_comm, &kernel_fn);
  tree->SetupFMM(&matrices);
  //
  pvfmm::PtFMM_Evaluate(tree, trg_value, n_trg);
  //
  delete tree;
}

// ---- custom SCTL stress kernels shared by the Stokes tractions S' and D' --------------------
// Stokes3D_DxT: double-layer source (density mu + source normal n) -> Cauchy stress tensor (9).
//   sigma_lm = (3/4pi)[ (2/3)(n.mu) d_lm/r^3 + (n_l r_p + a d_lp) r_m/r^5
//                       + (n_m r_p + a d_mp) r_l/r^5 - 10 a r_p r_l r_m/r^7 ] mu_p,  a = r.n.
// Verified term-by-term against the naive Sto3dDLPmat An (== fmm3d, 2.6e-15) and coded-vs-naive
// at 1.7e-14. u[p][l*3+m] holds the mu_p coefficient; ScaleFactor carries 3/(4pi). (Used by D'.)
struct pvfmmutils_Stokes3D_DxT_ {
  static const std::string& Name() { static const std::string s = "Stokes3D-DxT"; return s; }
  static constexpr sctl::Integer FLOPS() { return 60; }
  template <class Real> static constexpr Real uKerScaleFactor() { return 3 / (4 * sctl::const_pi<Real>()); }
  template <sctl::Integer digits, class VecType>
  static void uKerMatrix(VecType (&u)[3][9], const VecType (&r)[3], const VecType (&n)[3], const void*) {
    VecType r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    VecType rinv = sctl::approx_rsqrt<digits>(r2, r2 > VecType::Zero());
    VecType rinv2 = rinv*rinv, rinv3 = rinv2*rinv, rinv5 = rinv3*rinv2, rinv7 = rinv5*rinv2;
    const VecType one = (typename VecType::ScalarType)(1.0);
    const VecType twothird = (typename VecType::ScalarType)(2.0/3.0);
    const VecType ten = (typename VecType::ScalarType)(10.0);
    // DUMMY-NORMAL GUARD: SCTL wraps double-layer kernels as a pair; the "main" variant
    // (pvfmm's ker_poten, used ONLY by pvfmm's numeric scale/symmetry detection at kernel
    // Initialize) calls this with n = 0, which makes DxT identically ZERO. pvfmm's symmetry
    // search on a zero 3x9 kernel degenerates into a full 9^9 odometer enumeration (~160s PER
    // kernel object, ~8 min per DLPn call, single-threaded; profiled InitKernel=482s).
    // Substitute the RADIAL unit vector for the normal in that case: n_eff = n + (1-n.n) rhat.
    //   detection (n=0):  n_eff = rhat -> DxT(r, rhat): nonzero, r^-3 homogeneous, SAME
    //                     transformation law as the true DxT => correct scal/perm metadata,
    //                     detection prunes in ~ms. (A wrong-homogeneity stand-in, e.g. an
    //                     FxT-like r^-2 term, corrupts the far field to ~3e-2 -- the top
    //                     kernel's detected scaling IS used by the level machinery.)
    //   runtime (|n|=1):  n_eff = n up to 1e-16 rounding -> physics unchanged.
    VecType dummy = one - (n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);   // 1 for n=0 (detection), ~0 for unit n
    VecType ne0 = n[0] + dummy*r[0]*rinv;
    VecType ne1 = n[1] + dummy*r[1]*rinv;
    VecType ne2 = n[2] + dummy*r[2]*rinv;
    VecType ne[3] = {ne0, ne1, ne2};
    VecType a = r[0]*ne0 + r[1]*ne1 + r[2]*ne2;
    for (sctl::Integer p = 0; p < 3; p++)
      for (sctl::Integer l = 0; l < 3; l++)
        for (sctl::Integer m = 0; m < 3; m++) {
          VecType dlm = (typename VecType::ScalarType)(l==m ? 1.0 : 0.0);
          VecType dlp = (typename VecType::ScalarType)(l==p ? 1.0 : 0.0);
          VecType dmp = (typename VecType::ScalarType)(m==p ? 1.0 : 0.0);
          u[p][l*3+m] = twothird * ne[p] * dlm * rinv3
                      + (ne[l]*r[p] + a*dlp) * r[m] * rinv5
                      + (ne[m]*r[p] + a*dmp) * r[l] * rinv5
                      - ten * a * r[p] * r[l] * r[m] * rinv7;
        }
  }
};
using pvfmmutils_Stokes3D_DxT = sctl::GenericKernel<pvfmmutils_Stokes3D_DxT_>;

// Stokes3D_FSxT: [force(3), source/sink(1)] -> Cauchy stress (9). Force part = Stokeslet stress
// (== FxT, -6 r_i r_j r_k/r^5); source/sink part: the source flow u_i = C r_i/r^3 is div-free &
// harmonic => pressure 0 => stress = 2(d_jk/r^3 - 3 r_j r_k/r^5). Scale 1/(8pi), matching FSxU.
// This is the STRESS analogue of FSxU -- required so the M2M/M2L/L2L multipole is STRESS-consistent
// (its check potential is stress). Reading stress off a velocity multipole is lossy (~1e-2); a
// stress-consistent multipole gives ~1e-9 (verified: SLPn 1.3e-9, DLPn 4.0e-10 vs naive at acc 10).
struct pvfmmutils_Stokes3D_FSxT_ {
  static const std::string& Name() { static const std::string s = "Stokes3D-FSxT"; return s; }
  static constexpr sctl::Integer FLOPS() { return 50; }
  template <class Real> static constexpr Real uKerScaleFactor() { return 1 / (8 * sctl::const_pi<Real>()); }
  template <sctl::Integer digits, class VecType>
  static void uKerMatrix(VecType (&u)[4][9], const VecType (&r)[3], const void*) {
    VecType r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    VecType rinv = sctl::approx_rsqrt<digits>(r2, r2 > VecType::Zero());
    VecType rinv2 = rinv*rinv, rinv3 = rinv2*rinv, rinv5 = rinv3*rinv2;
    const VecType six = (typename VecType::ScalarType)(6.0);
    const VecType two = (typename VecType::ScalarType)(2.0);
    const VecType three = (typename VecType::ScalarType)(3.0);
    for (sctl::Integer i = 0; i < 3; i++)
      for (sctl::Integer j = 0; j < 3; j++)
        for (sctl::Integer k = 0; k < 3; k++)
          u[i][j*3+k] = -six * r[i]*r[j]*r[k]*rinv5;
    for (sctl::Integer j = 0; j < 3; j++)
      for (sctl::Integer k = 0; k < 3; k++) {
        VecType djk = (typename VecType::ScalarType)(j==k ? 1.0 : 0.0);
        u[3][j*3+k] = two*(djk*rinv3 - three*r[j]*r[k]*rinv5);
      }
  }
};
using pvfmmutils_Stokes3D_FSxT = sctl::GenericKernel<pvfmmutils_Stokes3D_FSxT_>;

// Stokes SLP velocity via SCTL ParticleFMM (unified with DLP: same velocity multipole, FSxU/FxU).
// NOTE: this replaces the raw-PtFMM StokesKernel::velocity() path, which was ~2x slower because
// velocity() carries a double-layer companion (sym_dip) whose translation matrices get built too;
// the pure single-layer FxU here builds a lighter operator (Precomp_Stokes3D-FxU_m10.data).
void stokes_slp_fmm( const vec&  sl_coord, const vec& sl_den,
                     const vec& trg_coord, vec& trg_value,
                     const int mult_order, const size_t max_pts, const size_t mem_size,
                     MPI_Comm& comm){
  (void)max_pts; (void)mem_size;
  const size_t n_sl  = sl_coord.size()  / PVFMM_COORD_DIM;
  const size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  sctl::Comm sctl_comm(comm);

  sctl::Stokes3D_FxU  kernel_sl;    // single-layer Stokeslet -> velocity
  sctl::Stokes3D_FSxU kernel_m2l;   // scale-invariant translations

  sctl::ParticleFMM<double,3> fmm(sctl_comm);
  fmm.SetAccuracy(static_cast<sctl::Integer>(mult_order));
  fmm.SetKernels(kernel_m2l, kernel_m2l, kernel_sl);   // (M2M, M2L, L2L)
  fmm.AddSrc("SL",  kernel_sl,  kernel_sl);             // (S2M, S2L)
  fmm.AddTrg("Vel", kernel_m2l, kernel_sl);             // (M2T, L2T)
  fmm.SetKernelS2T("SL", "Vel", kernel_sl);             // near

  sctl::Vector<double> src_coord(n_sl*3), src_den(n_sl*3), trg_c(n_trg*3);
  for (size_t i = 0; i < n_sl*3;  ++i) { src_coord[i] = sl_coord[i]; src_den[i] = sl_den[i]; }
  for (size_t i = 0; i < n_trg*3; ++i)   trg_c[i] = trg_coord[i];

  fmm.SetSrcCoord("SL", src_coord);
  fmm.SetSrcDensity("SL", src_den);
  fmm.SetTrgCoord("Vel", trg_c);

  sctl::Vector<double> U;
  fmm.Eval(U, "Vel");

  trg_value.assign(n_trg*3, 0.0);
  const size_t nout = std::min(static_cast<size_t>(U.Dim()), n_trg*3);
  for (size_t i = 0; i < nout; ++i) trg_value[i] = U[i];
}

// Stokes SLP traction S' (adjoint DLP): t_i = sigma_ij n^t_j, sigma = Cauchy stress of the
// single-layer (Stokeslet) velocity field. Via SCTL ParticleFMM with a STRESS-CONSISTENT wiring
// (unified with the DLP traction D'): single-layer source -> stress (FxT), stress translations
// (FSxT/FxT) so the multipole reproduces the stress field, then contract the 9-component stress
// with the target normal. S2T = FxT => uses Precomp_Stokes3D-FxT_m10.data.
void stokes_slpn_fmm( const vec&  sl_coord, const vec& sl_den,
                      const vec& trg_coord, const vec& trg_normal, vec& trg_value,
                      const int mult_order, const size_t max_pts, const size_t mem_size,
                      MPI_Comm& comm){
  (void)max_pts; (void)mem_size;
  const size_t n_sl  = sl_coord.size()  / PVFMM_COORD_DIM;
  const size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  sctl::Comm sctl_comm(comm);

  sctl::Stokes3D_FxT       kernel_ft;    // Stokeslet -> stress (S2M/S2L/S2T/L2L/L2T)
  pvfmmutils_Stokes3D_FSxT kernel_fst;   // [force+source/sink] -> stress (M2M/M2L/M2T)

  sctl::ParticleFMM<double,3> fmm(sctl_comm);
  fmm.SetAccuracy(static_cast<sctl::Integer>(mult_order));
  fmm.SetKernels(kernel_fst, kernel_fst, kernel_ft);   // (M2M, M2L, L2L) -> stress-consistent multipole
  fmm.AddSrc("SL", kernel_ft, kernel_ft);              // (S2M, S2L) = Stokeslet -> stress
  fmm.AddTrg("Stress", kernel_fst, kernel_ft);         // (M2T, L2T) -> stress
  fmm.SetKernelS2T("SL", "Stress", kernel_ft);         // near = Stokeslet -> stress

  sctl::Vector<double> src_coord(n_sl*3), src_den(n_sl*3), trg_c(n_trg*3);
  for (size_t i = 0; i < n_sl*3;  ++i) { src_coord[i] = sl_coord[i]; src_den[i] = sl_den[i]; }
  for (size_t i = 0; i < n_trg*3; ++i)   trg_c[i] = trg_coord[i];

  fmm.SetSrcCoord("SL", src_coord);
  fmm.SetSrcDensity("SL", src_den);
  fmm.SetTrgCoord("Stress", trg_c);

  sctl::Vector<double> S;
  fmm.Eval(S, "Stress");                          // 9-component stress per target

  trg_value.assign(n_trg*3, 0.0);
  const size_t navail = static_cast<size_t>(S.Dim()) / 9;
  for (size_t q = 0; q < n_trg && q < navail; ++q) {
    const double* nt = &trg_normal[q*3];
    for (int l = 0; l < 3; ++l)
      trg_value[q*3+l] = S[q*9+l*3+0]*nt[0] + S[q*9+l*3+1]*nt[1] + S[q*9+l*3+2]*nt[2];
  }
}

// Stokes DLP (double-layer / stresslet velocity) via SCTL's ParticleFMM (pvfmm-backed).
//
// pvfmm's raw PtFMM double-layer (dl_coord/dl_den) mechanism does NOT correctly evaluate the
// Stokes double-layer: it auto-generates all translation operators from one kernel, and that
// auto-generation is inconsistent for the (non-scale-invariant) double-layer far field
// (all-direct vs multilevel self-consistency relerr ~ O(1); verified in scratchpad). The
// author's own SCTL/CSBQ stack does it differently and CORRECTLY: ParticleFMM lets each
// translation operator use a DIFFERENT kernel -- the source->multipole / source->target use
// the double-layer kernel (Stokes3D_DxU), while the M2M/M2L/L2L multipole machinery uses the
// scale-invariant Stokes3D_FSxU. This mirrors ParticleFMM<Real,3>::test() (fmm-wrapper.txx).
//
// Stokes3D_DxU is exactly the standard stresslet u_i = (3/4pi) r_i r_j (r.n)/r^5 * mu_j,
// matching fmm3d's stfmm3d stresslet (T_ijk = 3 r_i r_j r_k/(4pi r^5)); density mu and the
// source normal n are supplied separately (SetSrcCoord carries the normal). Validated in the
// test against fmm3d's Sto3dDLPfmm_il and against ParticleFMM::EvalDirect.
//   dl_coord  : source coords, node-interleaved [x1 y1 z1 x2 ...]  (3*n_dl)
//   dl_sigma  : DLP density mu, node-interleaved                    (3*n_dl)
//   dl_normal : source normals nu, node-interleaved                 (3*n_dl)
//   mult_order: interpreted as SCTL accuracy DIGITS (SetAccuracy), not a pvfmm multipole order
void stokes_dlp_fmm( const vec&  dl_coord, const vec& dl_sigma, const vec& dl_normal,
                     const vec& trg_coord, vec& trg_value,
                     const int mult_order, const size_t max_pts, const size_t mem_size,
                     MPI_Comm& comm){
  (void)max_pts; (void)mem_size;   // ParticleFMM manages its own tree via SetAccuracy
  const size_t n_dl  = dl_coord.size()  / PVFMM_COORD_DIM;
  const size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  sctl::Comm sctl_comm(comm);

  sctl::Stokes3D_FxU  kernel_sl;    // single-layer Stokeslet (F x U)
  sctl::Stokes3D_DxU  kernel_dl;    // double-layer stresslet (D x U): (3/4pi) r_i r_j (r.n)/r^5
  sctl::Stokes3D_FSxU kernel_m2l;   // scale-invariant single-layer+source/sink -> velocity,
                                    // required for the M2M/M2L/L2L/M2T translations

  sctl::ParticleFMM<double,3> fmm(sctl_comm);
  fmm.SetAccuracy(static_cast<sctl::Integer>(mult_order));
  fmm.SetKernels(kernel_m2l, kernel_m2l, kernel_sl);   // (M2M, M2L, L2L)
  fmm.AddSrc("DL",  kernel_dl,  kernel_dl);             // (S2M, S2L)
  fmm.AddTrg("Vel", kernel_m2l, kernel_sl);             // (M2T, L2T)
  fmm.SetKernelS2T("DL", "Vel", kernel_dl);             // near / source->target

  sctl::Vector<double> src_coord(n_dl*3), src_normal(n_dl*3), src_den(n_dl*3), trg_c(n_trg*3);
  for (size_t i = 0; i < n_dl*3;  ++i) { src_coord[i] = dl_coord[i]; src_normal[i] = dl_normal[i]; src_den[i] = dl_sigma[i]; }
  for (size_t i = 0; i < n_trg*3; ++i)   trg_c[i] = trg_coord[i];

  fmm.SetSrcCoord("DL", src_coord, src_normal);   // source coords + normals
  fmm.SetSrcDensity("DL", src_den);               // DLP density mu (3 per source)
  fmm.SetTrgCoord("Vel", trg_c);

  sctl::Vector<double> U;
  fmm.Eval(U, "Vel");                             // pvfmm-backed far field

  trg_value.assign(n_trg*3, 0.0);
  const size_t nout = std::min(static_cast<size_t>(U.Dim()), n_trg*3);
  for (size_t i = 0; i < nout; ++i) trg_value[i] = U[i];
}

// (Stokes3D_DxT and Stokes3D_FSxT kernel structs are defined above stokes_slpn_fmm -- shared by
//  the single-layer traction S' and the double-layer traction D'.)

// Stokes DLP TRACTION (D') via SCTL ParticleFMM, STRESS-CONSISTENT wiring: the M2M/M2L/L2L
// translations produce STRESS check potentials (FSxT/FxT) so the multipole reproduces the stress
// field (mirrors pvfmm PtFMM StokesKernel::stress), NOT a velocity multipole read as stress.
//   S2M/S2L = DxT (double-layer -> stress), M2M/M2L = FSxT, L2L = FxT, M2T = FSxT, L2T = FxT,
//   S2T = DxT.  Output = 9-comp stress per target; contract with the TARGET normal -> traction.
// Verified vs naive stress-of-double-layer (== fmm3d Sto3dDLPnfmm_il): 4.0e-10 at acc 10, N>=40000.
//   dl_coord/dl_sigma/dl_normal : source coords / DLP density mu / source normals (3*n_dl each)
//   trg_normal : target UNIT normals (3*n_trg), contracted with the stress
//   mult_order : SCTL accuracy DIGITS
void stokes_dlpn_fmm( const vec&  dl_coord, const vec& dl_sigma, const vec& dl_normal,
                      const vec& trg_coord, const vec& trg_normal, vec& trg_value,
                      const int mult_order, const size_t max_pts, const size_t mem_size,
                      MPI_Comm& comm){
  (void)max_pts; (void)mem_size;
  const size_t n_dl  = dl_coord.size()  / PVFMM_COORD_DIM;
  const size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  sctl::Comm sctl_comm(comm);

  sctl::Stokes3D_FxT       kernel_ft;    // Stokeslet -> stress (L2L, L2T)
  pvfmmutils_Stokes3D_FSxT kernel_fst;   // [force+source/sink] -> stress (M2M, M2L, M2T)
  pvfmmutils_Stokes3D_DxT  kernel_dt;    // double-layer -> stress (S2M, S2L, S2T)

  sctl::ParticleFMM<double,3> fmm(sctl_comm);
  fmm.SetAccuracy(static_cast<sctl::Integer>(mult_order));
  fmm.SetKernels(kernel_fst, kernel_fst, kernel_ft);   // (M2M, M2L, L2L) -> stress-consistent multipole
  fmm.AddSrc("DL", kernel_dt, kernel_dt);              // (S2M, S2L) = double-layer -> stress
  fmm.AddTrg("Stress", kernel_fst, kernel_ft);         // (M2T, L2T) -> stress
  fmm.SetKernelS2T("DL", "Stress", kernel_dt);         // near = double-layer -> stress

  sctl::Vector<double> src_coord(n_dl*3), src_normal(n_dl*3), src_den(n_dl*3), trg_c(n_trg*3);
  for (size_t i = 0; i < n_dl*3;  ++i) { src_coord[i] = dl_coord[i]; src_normal[i] = dl_normal[i]; src_den[i] = dl_sigma[i]; }
  for (size_t i = 0; i < n_trg*3; ++i)   trg_c[i] = trg_coord[i];

  fmm.SetSrcCoord("DL", src_coord, src_normal);
  fmm.SetSrcDensity("DL", src_den);
  fmm.SetTrgCoord("Stress", trg_c);

  sctl::Vector<double> S;
  fmm.Eval(S, "Stress");                          // 9-component stress per target

  // contract stress with the target normal: t_l = sigma_lm n^t_m
  trg_value.assign(n_trg*3, 0.0);
  const size_t navail = static_cast<size_t>(S.Dim()) / 9;
  for (size_t q = 0; q < n_trg && q < navail; ++q) {
    const double* nt = &trg_normal[q*3];
    for (int l = 0; l < 3; ++l)
      trg_value[q*3+l] = S[q*9+l*3+0]*nt[0] + S[q*9+l*3+1]*nt[1] + S[q*9+l*3+2]*nt[2];
  }
}

// ---- custom SCTL kernel for the Laplace DLP normal-derivative (D', hypersingular) -----------
// Laplace3D_DxdU: dipole source (charge sigma + source normal n) -> gradient (3). The DLP potential
// is (1/4pi)(r.n)/r^3 * sigma; its gradient is (1/4pi)[ n_i/r^3 - 3 r_i (r.n)/r^5 ] sigma.
// KDIM0=1 (charge), KDIM1=3 (gradient), normal-dim=3. Verified coded-vs-naive 2.3e-14.
struct pvfmmutils_Laplace3D_DxdU_ {
  static const std::string& Name() { static const std::string s = "Laplace3D-DxdU"; return s; }
  static constexpr sctl::Integer FLOPS() { return 20; }
  template <class Real> static constexpr Real uKerScaleFactor() { return 1 / (4 * sctl::const_pi<Real>()); }
  template <sctl::Integer digits, class VecType>
  static void uKerMatrix(VecType (&u)[1][3], const VecType (&r)[3], const VecType (&n)[3], const void*) {
    VecType r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    VecType rinv = sctl::approx_rsqrt<digits>(r2, r2 > VecType::Zero());
    VecType rinv3 = rinv*rinv*rinv, rinv5 = rinv3*rinv*rinv;
    VecType rn = r[0]*n[0]+r[1]*n[1]+r[2]*n[2];
    const VecType three = (typename VecType::ScalarType)(3.0);
    for (sctl::Integer i = 0; i < 3; i++) u[0][i] = n[i]*rinv3 - three*r[i]*rn*rinv5;
  }
};
using pvfmmutils_Laplace3D_DxdU = sctl::GenericKernel<pvfmmutils_Laplace3D_DxdU_>;

// Laplace DLPn (normal-derivative of the double-layer) via SCTL ParticleFMM, GRADIENT-CONSISTENT
// wiring: the translations produce GRADIENT check potentials (Laplace3D_FxdU) so the multipole
// reproduces the gradient field, then the 3-component gradient is contracted with the TARGET normal.
//   S2M/S2L = DxdU, M2M/M2L/L2L = FxdU, M2T/L2T = FxdU, S2T = DxdU.
// (Laplace needs no source/sink augmentation -- unlike Stokes.) Verified vs naive: 5.0e-10 at acc 10.
//   dl_coord/dl_normal : source coords / source normals (3*n_dl);  dl_sigma : charge (n_dl)
//   trg_normal : target UNIT normals (3*n_trg);  trg_value : DLPn scalar per target (n_trg)
void laplace_dlpn_fmm( const vec&  dl_coord, const vec& dl_normal, const vec& dl_sigma,
                       const vec& trg_coord, const vec& trg_normal, vec& trg_value,
                       const int mult_order, const size_t max_pts, const size_t mem_size,
                       MPI_Comm& comm){
  (void)max_pts; (void)mem_size;
  const size_t n_dl  = dl_coord.size()  / PVFMM_COORD_DIM;
  const size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  sctl::Comm sctl_comm(comm);

  sctl::Laplace3D_FxdU     kernel_fd;    // charge -> gradient (translations, L2T)
  pvfmmutils_Laplace3D_DxdU kernel_dd;   // dipole -> gradient (S2M/S2L/S2T, custom)

  sctl::ParticleFMM<double,3> fmm(sctl_comm);
  fmm.SetAccuracy(static_cast<sctl::Integer>(mult_order));
  fmm.SetKernels(kernel_fd, kernel_fd, kernel_fd);   // (M2M,M2L,L2L) -> gradient-consistent multipole
  fmm.AddSrc("DL", kernel_dd, kernel_dd);            // (S2M,S2L) = dipole -> gradient
  fmm.AddTrg("Grad", kernel_fd, kernel_fd);          // (M2T,L2T) -> gradient
  fmm.SetKernelS2T("DL", "Grad", kernel_dd);         // near = dipole -> gradient

  sctl::Vector<double> src_coord(n_dl*3), src_normal(n_dl*3), src_den(n_dl), trg_c(n_trg*3);
  for (size_t i = 0; i < n_dl*3;  ++i) { src_coord[i] = dl_coord[i]; src_normal[i] = dl_normal[i]; }
  for (size_t i = 0; i < n_dl;    ++i)   src_den[i] = dl_sigma[i];
  for (size_t i = 0; i < n_trg*3; ++i)   trg_c[i] = trg_coord[i];

  fmm.SetSrcCoord("DL", src_coord, src_normal);
  fmm.SetSrcDensity("DL", src_den);
  fmm.SetTrgCoord("Grad", trg_c);

  sctl::Vector<double> G;
  fmm.Eval(G, "Grad");                            // 3-component gradient per target

  trg_value.assign(n_trg, 0.0);
  const size_t navail = static_cast<size_t>(G.Dim()) / 3;
  for (size_t q = 0; q < n_trg && q < navail; ++q) {
    const double* nt = &trg_normal[q*3];
    trg_value[q] = G[q*3+0]*nt[0] + G[q*3+1]*nt[1] + G[q*3+2]*nt[2];
  }
}

//
void laplace_slpn_fmm( const vec&  sl_coord, const vec& sl_den,
                       const vec& trg_coord, const vec& trg_normal, vec& trg_value,
                       const int mult_order, const size_t max_pts, const size_t mem_size,
                       MPI_Comm& comm){
  const pvfmm::Kernel<double>& kernel_fn = pvfmm::LaplaceKernel<double>::gradient();
  const size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  //
  sctl::Comm sctl_comm(comm);
  // size_t max_pts = 600;
  // Empty double-layer inputs for SLP-only run.
  vec dl_coord, dl_den;
  vec trg_grad_value(n_trg * 3);
  auto* tree = PtFMM_CreateTree(sl_coord, sl_den, dl_coord, dl_den,
                                trg_coord, sctl_comm, (int)max_pts, pvfmm::FreeSpace);
  //
  pvfmm::PtFMM<double> matrices;
  matrices.Initialize(mult_order, sctl_comm, &kernel_fn);
  tree->SetupFMM(&matrices);
  //
  pvfmm::PtFMM_Evaluate(tree, trg_grad_value, n_trg);
  //
  delete tree;
  //
  trg_value.assign(n_trg, 0.0);
  for (size_t i = 0; i < n_trg; ++i) {
    const size_t gi = i * PVFMM_COORD_DIM;
    trg_value[i] = trg_grad_value[gi + 0] * trg_normal[gi + 0]
                 + trg_grad_value[gi + 1] * trg_normal[gi + 1]
                 + trg_grad_value[gi + 2] * trg_normal[gi + 2];
  }
}

// kernel_fn.dbl_layer_poten --> DLP in kernel.txx
void laplace_dlp_fmm( const vec&  dl_coord, const vec& dl_normal, const vec& dl_sigma,
                      const vec& trg_coord, vec& trg_value,
                      const int mult_order, const size_t max_pts, const size_t mem_size,
                      MPI_Comm& comm){
  const pvfmm::Kernel<double>& kernel_fn = pvfmm::LaplaceKernel<double>::potential();
  const size_t n_src = dl_coord.size() / PVFMM_COORD_DIM;
  const size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  //
  sctl::Comm sctl_comm(comm);
  // size_t max_pts = 600;
  // Prepare DLP density as [nx, ny, nz, sigma] per source point.
  vec dl_den(n_src * 4);
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n_src; ++i) {
    dl_den[4*i + 0] = dl_normal[3*i + 0];
    dl_den[4*i + 1] = dl_normal[3*i + 1];
    dl_den[4*i + 2] = dl_normal[3*i + 2];
    dl_den[4*i + 3] = dl_sigma[i];
  }
  // Empty single-layer inputs for DLP-only run.
  vec sl_coord, sl_den;
  auto* tree = PtFMM_CreateTree(sl_coord, sl_den, dl_coord, dl_den,
                                trg_coord, sctl_comm, (int)max_pts, pvfmm::FreeSpace);
  //
  pvfmm::PtFMM<double> matrices;
  matrices.Initialize(mult_order, sctl_comm, &kernel_fn);
  tree->SetupFMM(&matrices);
  //
  pvfmm::PtFMM_Evaluate(tree, trg_value, n_trg);
  //
  delete tree;
}

// SLP (fixed k in helmholtz_poten_).
void helmholtz_slp_fmm( const vec&  sl_coord, const vec& sl_den,
                        const vec& trg_coord, vec& trg_value,
                        const int mult_order, const size_t max_pts, const size_t mem_size,
                        MPI_Comm& comm){
  const pvfmm::Kernel<double>& kernel_fn = pvfmm::HelmholtzKernel<double>::potential();
  const size_t n_trg = trg_coord.size() / PVFMM_COORD_DIM;
  //
  sctl::Comm sctl_comm(comm);
  // Empty double-layer inputs for SLP-only run.
  vec dl_coord, dl_den;
  auto* tree = PtFMM_CreateTree(sl_coord, sl_den, dl_coord, dl_den,
                                trg_coord, sctl_comm, (int)max_pts, pvfmm::FreeSpace);
  //
  pvfmm::PtFMM<double> matrices;
  matrices.Initialize(mult_order, sctl_comm, &kernel_fn);
  tree->SetupFMM(&matrices);
  //
  pvfmm::PtFMM_Evaluate(tree, trg_value, n_trg);
  //
  delete tree;
}



#endif // PVFMM_WRAP_H
