# Finite Difference Method
### Main Scripts
* `elliptic.m:` Solves 2D elliptic PDE using 5-point and 9-point molecule 
* `parabolic_fw_bw.m:` Solves 1D parabolic PDE using forward or backward scheme
* `parabolic_cn.m:` Solves 1D parabolic PDE using Crank-Nicolson method
* `hyperbolic_explicit.m:` Solves 1D hyperbolic PDE using explicit scheme
* `hyperbolic_implicit.m:` Solves 1D hyperbolic PDE using implicit scheme
### Utility Functions
* `five_point.m:` Creates coefficient matrix representing the 5-point stencil (2D case)
* `nine_point.m:` Creates coefficient matrix representing the 9-point stencil (2D case)
* `poisson3D.m:` Creates coefficient matrix representing the 7-point stencil (3D case)
* `convdiff.m:` Creates a simple convection-diffusion matrix for the 2D or 3D case
* `convdiff2D.m:` Creates a more complex convection-diffusion matrix for the 2D case
* `convdiff3D.m:` Creates a more complex convection-diffusion matrix for the 3D case
