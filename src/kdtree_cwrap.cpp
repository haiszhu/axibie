// kdtree_cwrap.cpp
//
// Single-file, multi-backend extern "C" bridge: build nbodyhpc and taiya kd-trees and
// call them from Fortran (kdtree_mod), so the two libraries can be compared.
//
//   backend "nbody" -> haiszhu/nbodyhpc (external/nbodyhpc): 3D, float32, AVX2/NEON.
//                      No native ball query -> BUILD + EXPORT here; radius traversal in
//                      Fortran (kdtree_mod::ball_bbox_query).
//   backend "taiya" -> taiya/kdtree (external/kdtree): double, N-dim, native C++ ball.
//                      BUILD + QUERY + return indices, all here.
//
// NAME-CLASH: both define a class `KDTree`. nbodyhpc -> wenda::kdtree::KDTree (we always
// fully-qualify and never `using` it); taiya -> ::KDTree. They never collide.
//
// TAIYA HEADERS ARE BUILD-PATCHED (see Makefile_kdtree), submodule stays pristine:
// MyHeaps.h has its C++17-removed `throw(HeapEmptyException)` specs stripped so the whole
// file builds at C++20. KDTree.h is used as-is: it calls mexErrMsgTxt/mexPrintf outside
// any #ifdef, so this file is compiled with the MATLAB SDK include and linked with libmex.
//
// Boundary (FMM3D/vec-kernels, -i8): symbols end in "_", args by pointer; int64_t <->
// integer*8, float <-> real*4, double <-> real*8.

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

#include <kdtree/kdtree.hpp> // nbodyhpc -> wenda::kdtree::KDTree  (include FIRST: C++20, clean)
#include "KDTree.h"          // taiya    -> ::KDTree  (patched copy from $(GEN); `using namespace std`)

// ============================================================================
//  Backend: nbodyhpc  (3D float32; build + export, ball traversal in Fortran)
// ============================================================================

namespace {
thread_local std::unique_ptr<wenda::kdtree::KDTree> g_nbody;
} // namespace

extern "C" {

void kdtree_nbody_build_(const int64_t *np, const float *pts,
                         const int64_t *leafsize, const int64_t *maxthreads,
                         int64_t *nnodes, int64_t *npad) {
    auto const *p = reinterpret_cast<std::array<float, 3> const *>(pts);
    tcb::span<const std::array<float, 3>> positions(p, static_cast<size_t>(*np));

    wenda::kdtree::KDTreeConfiguration cfg;
    cfg.leaf_size = static_cast<int>(*leafsize);
    cfg.max_threads = static_cast<int>(*maxthreads);

    g_nbody = std::make_unique<wenda::kdtree::KDTree>(positions, cfg);

    *nnodes = static_cast<int64_t>(g_nbody->nodes().size());
    *npad = static_cast<int64_t>(g_nbody->positions().size());
}

// Export the retained tree into Fortran-allocated buffers (see kdtree_mod for layout).
void kdtree_nbody_export_(int64_t *dim, float *split, int64_t *lft, int64_t *rgt,
                          float *xp, float *yp, float *zp, int64_t *idx) {
    auto nodes = g_nbody->nodes();
    for (size_t i = 0; i < nodes.size(); ++i) {
        dim[i] = nodes[i].dimension_;
        split[i] = nodes[i].split_;
        lft[i] = static_cast<int64_t>(nodes[i].left_);
        rgt[i] = static_cast<int64_t>(nodes[i].right_);
    }
    auto const &pos = g_nbody->positions();
    size_t n = pos.size();
    for (size_t s = 0; s < n; ++s) {
        xp[s] = pos.positions_[0][s];
        yp[s] = pos.positions_[1][s];
        zp[s] = pos.positions_[2][s];
        idx[s] = static_cast<int64_t>(pos.indices_[s]);
    }
}

void kdtree_nbody_free_() { g_nbody.reset(); }

} // extern "C"

// ============================================================================
//  Backend: taiya/kdtree  (double, N-dim; native C++ ball query)
// ============================================================================

namespace {
thread_local std::unique_ptr<KDTree> g_taiya; // ::KDTree (taiya)
} // namespace

extern "C" {

// Build a taiya tree over `np` points (interleaved xyz, length 3*np, double).
void kdtree_taiya_build_(const int64_t *np, const double *pts) {
    std::vector<Point> P(static_cast<size_t>(*np), Point(3));
    for (int64_t i = 0; i < *np; ++i) {
        P[i][0] = pts[3 * i + 0];
        P[i][1] = pts[3 * i + 1];
        P[i][2] = pts[3 * i + 2];
    }
    g_taiya = std::make_unique<KDTree>(P);
}

// Ball (radius) query on the retained tree: original 0-based indices within `radius` of
// `qpt` (length 3); up to `nmax` written to `out`, count in `nout`. Uses taiya's ball_query.
void kdtree_taiya_ball_(const double *qpt, const double *radius, const int64_t *nmax,
                        int64_t *out, int64_t *nout) {
    Point q(3);
    q[0] = qpt[0];
    q[1] = qpt[1];
    q[2] = qpt[2];

    std::vector<int> idxs;
    std::vector<double> dists;
    g_taiya->ball_query(q, *radius, idxs, dists);

    int64_t n = static_cast<int64_t>(idxs.size());
    int64_t nwrite = (n < *nmax) ? n : *nmax;
    for (int64_t i = 0; i < nwrite; ++i)
        out[i] = static_cast<int64_t>(idxs[i]); // 0-based original index
    *nout = n;                                   // true count (caller checks n <= nmax)
}

void kdtree_taiya_free_() { g_taiya.reset(); }

} // extern "C"
