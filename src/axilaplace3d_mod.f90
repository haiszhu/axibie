module axilaplace3d_mod
  ! ------------------------------------------------------------------
  ! Base module for the axisymmetric Laplace subsystem.
  !
  ! Thin RE-EXPORT WRAPPER: it brings the shared kind parameters
  ! (r64, r128, c64, c128) and the Gauss-Legendre primitive
  ! (gauss_r64, gauss_r128) from axistokes3d_mod and re-exports them under
  ! the Laplace name, so the Laplace modules depend on a Laplace-named base
  ! rather than the Stokes one:
  !   use axilaplace3d_mod, only: r64, ...
  ! Single source of truth -- the actual definitions live in axistokes3d_mod;
  ! this module adds no new definitions, only the Laplace-named handle.
  ! ------------------------------------------------------------------
  use axistokes3d_mod, only: r64, r128, c64, c128, gauss_r64, gauss_r128
  implicit none
  public :: r64, r128, c64, c128, gauss_r64, gauss_r128
end module axilaplace3d_mod
