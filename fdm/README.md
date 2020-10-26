# Finite Difference Method
### Main Scripts
* `elliptic.m:` Solves 2D elliptic PDE using 5-point and 9-point molecule 
* `parabolic_fw_bw.m:` Solves 1D parabolic PDE using forward or backward scheme
* `parabolic_cn.m:` Solves 1D parabolic PDE using Crank-Nicolson method
* `convdiff2d_stencil.m:` Solves 2D convection-diffusion PDE
* `hyperbolic_explicit.m:` Solves 1D hyperbolic PDE using explicit scheme
* `hyperbolic_implicit.m:` Solves 1D hyperbolic PDE using implicit scheme
### Utility Functions
* `five_point.m:` Creates coefficient matrix representing the 5-point stencil (2D case)
* `nine_point.m:` Creates coefficient matrix representing the 9-point stencil (2D case)
* `poisson3D.m:` Creates coefficient matrix representing the 7-point stencil (3D case)
* `convdiff.m:` Creates a simple convection-diffusion matrix for the 2D or 3D case
* `convdiff2D.m:` Creates a more complex convection-diffusion matrix for the 2D case (see reference [paper](https://hal.inria.fr/inria-00466828))
* `convdiff3D.m:` Creates a more complex convection-diffusion matrix for the 3D case (see reference [paper](https://hal.inria.fr/inria-00466828))

#### Reference
> Giraud, L., Haidar, A. and Saad, Y., 2010. Sparse approximations of the Schur complement for parallel algebraic hybrid linear solvers in 3D. 
[paper link](https://hal.inria.fr/inria-00466828)
