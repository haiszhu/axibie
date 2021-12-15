# axibie
Axisymmetric Boundary Integral Equation Tools in Matlab

The goal of this project is to provide tools in Matlab for solving the axisymmetric Laplace and Stokes problems. The numerical method in use is based on a high-order accurate panel-based boundary integral scheme. 

This project is in a very early state. It contains functions from Alex Barnett's package BIE2D available at: https://github.com/ahbarnett/BIE2D.

## Examples

## References

1. Guo, Hanliang, Hai Zhu, Ruowen Liu, Marc Bonnet, and Shravan Veerapaneni. 2021. “Optimal Slip Velocities of Micro-Swimmers with Arbitrary Axisymmetric Shapes.” *Journal of Fluid Mechanics* 910.

2. Hao, Sijia, Alex H Barnett, Per-Gunnar Martinsson, and P Young. 2014. “High-Order Accurate Methods for Nyström Discretization of Integral Equations on Smooth Curves in the Plane.” *Advances in Computational Mathematics* 40 (1): 245–72.

3. Helsing, Johan, and Anders Karlsson. 2014. “An Explicit Kernel-Split Panel-Based Nyström Scheme for Integral Equations on Axially Symmetric Surfaces.” *Journal of Computational Physics* 272: 686–703.


4. Veerapaneni, Shravan K, Denis Gueyffier, George Biros, and Denis Zorin. 2009. “A Numerical Method for Simulating the Dynamics of 3d Axisymmetric Vesicles Suspended in Viscous Flows.” *Journal of Computational Physics* 228 (19): 7233–49.

5. Wu, Bowei, Hai Zhu, Alex Barnett, and Shravan Veerapaneni. 2020. “Solution of Stokes Flow in Complex Nonsmooth 2d Geometries via a Linear-Scaling High-Order Adaptive Integral Equation Scheme.” *Journal of Computational Physics* 410: 109361.


## To do list

* Implement high order Fourier modes to enable nonsymmetric potential and flow simulation
* Implement multiple partciles + possible interaction with confined geometry
* Accelerate via FMM?

