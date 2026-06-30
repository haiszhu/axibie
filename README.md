# axibie
Axisymmetric Boundary Integral Equation Tools in MATLAB

This repository provides tools in Matlab for solving the axisymmetric Laplace and Stokes problems. The numerical method in use is based on a high-order accurate panel-based boundary integral scheme. 

This project is in an early state. It contains functions from [BIE2D](https://github.com/ahbarnett/BIE2D).

The code supports modal/toroidal Laplace and Stokes Green's function in [`/src`](src) in fortran. All numerical results in [`/test/stokes`](test/stokes) and [`/test/laplace`](test/laplace) have been verified for Stokes and Laplace boundary value problems. [`/test/stokes`](test/stokes) and [`/test/laplace`](test/laplace) contain MATLAB codes calling [`/src`](src) fortran using MEX gateway functions automatically created via [mwrap](https://github.com/zgimbutas/mwrap).

The quadrature for modal/toroidal Laplace and Stokes Green's function is high-order kernel-split (see references), allowing near machine precision accuracies for targets near/on the boundary via one unified interface. For part of the kernel-split formula, see references. For the implementation, it is still in progress.

## Installation

Prerequisites: MATLAB (with its MEX compiler SDK), a recent `gfortran` (the [`Makefile`](Makefile) uses `gfortran-15`), and [mwrap](https://github.com/zgimbutas/mwrap) (expected at `~/mwrap/mwrap`).

Build the Fortran compute layer and the MEX gateway from the repository root:

```sh
git submodule update --init --recursive
make -f Makefile_kdtree
make mex
```

You can now run any driver in [`/test/stokes`](test/stokes) or [`/test/laplace`](test/laplace), e.g.

```matlab
run('test/laplace/test_axissymslap_lap_slp_bvp_0th.m')
run('test/laplace/test_axissymslap_lap_slp_bvp.m')
```

## Examples

A minimal close-evaluation: build a panel curve, then assemble the single-layer block matrix for a target sitting just outside the surface (where naive quadrature loses digits) — accurate to near machine precision.

```matlab
addpath utils matlab                               % from the repo root, after `make mex`

p = 16; np = 14;                                   % 14 panels x 16 Gauss nodes each
s = []; s.p = p; s.Z = @(t) sin(t) - 1i*cos(t);    % a curve: unit-sphere meridian (rho = sin t, z = -cos t)
s.tpan = linspace(0,pi,np+1)'; s = quadr(s,[],'p','G');

d = 1e-3;                                          % a target a distance d just OUTSIDE the surface
t = []; t.x = (1+d)*(sin(1.2) - 1i*cos(1.2));      %   (a unit-sphere point pushed out along its normal)

% 0th-mode axisymmetric Laplace single-layer block (iside=1 exterior, iclosed=0 open arc with axis poles)
A = axls_slp_blockmat_mex(numel(t.x), t.x, p, np, ...
      s.x, s.nx, s.ws, s.wxp, s.tpan, s.xlo, s.xhi, 1, 0, []);

u = A * ones(numel(s.x),1);                        % single-layer potential of a unit density at t
```

`A` is the `nt x N` close-evaluation operator (here `1 x 224`); swap `axls_slp_blockmat_mex` for `axls_{slpn,dlp,dlpn}_blockmat_mex` (Laplace) or `axss_*` (Stokes), and the `_nmode_mex` variants for all azimuthal modes at once. See [`/test/stokes`](test/stokes) and [`/test/laplace`](test/laplace) for full boundary-value-problem drivers.

## References

1. The Pozrikidis

2. Hanliang Guo, Hai Zhu, Ruowen Liu, Marc Bonnet, and Shravan Veerapaneni. 2021. “Optimal Slip Velocities of Micro-Swimmers with Arbitrary Axisymmetric Shapes.” *Journal of Fluid Mechanics* 910.

3. Sijia Hao, Alex H Barnett, Per-Gunnar Martinsson, and P Young. 2014. “High-Order Accurate Methods for Nyström Discretization of Integral Equations on Smooth Curves in the Plane.” *Advances in Computational Mathematics* 40 (1): 245–72.

4. Johan Helsing, and Anders Karlsson. 2014. “An Explicit Kernel-Split Panel-Based Nyström Scheme for Integral Equations on Axially Symmetric Surfaces.” *Journal of Computational Physics* 272: 686–703.

5. Shravan Veerapaneni, Denis Gueyffier, George Biros, and Denis Zorin. 2009. “A Numerical Method for Simulating the Dynamics of 3d Axisymmetric Vesicles Suspended in Viscous Flows.” *Journal of Computational Physics* 228 (19): 7233–49.

6. Bowei Wu, Hai Zhu, Alex Barnett, and Shravan Veerapaneni. 2020. “Solution of Stokes Flow in Complex Nonsmooth 2d Geometries via a Linear-Scaling High-Order Adaptive Integral Equation Scheme.” *Journal of Computational Physics* 410: 109361.

7. James Garritano, Yuval Kluger, Vladimir Rokhlin, Kirill Serkh. 2022. "On the efficient evaluation of the azimuthal Fourier components of the Green's function for Helmholtz's equation in cylindrical coordinates." *Journal of Computational Physics* 471: 111585.


## To do list

* (debug DLPn done, need more testing) Implement high order Fourier modes to enable nonsymmetric potential and flow simulation
* (prototype done, need cleanup and memory footprint reduction) Implement multiple partciles + possible interaction with confined geometry
* (what is SLPnn, DLPnn? Azz in Dspecialquad for 2D?) Laplace case is a by-product, why not
* (prototype done) Accelerate via FMM3D/PVFMM/DMK
* (need to learn) 3d flow visualization
* (do we need variation in centerline/not axisymmetric any more) CSBQ
