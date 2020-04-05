# Finite Element Method
### Main Scripts
* `fem3.m:` Solves 2D Helmholtz equation on a square domain using triangular finite elements
* `fem4.m:` Solves 2D Helmholtz equation on a square domain using square finite elements
* `fem_helmholtz.m:` Solves 2D Helmholtz equation on a square domain using triangular or square finite elements, and benchmarks various Krylov solvers combined with Jacobi, Gauss-Seidel, and Symmetric Gauss-Seidel preconditioners (see chapter 2 in [doc](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf) for more details)
* `helm_feti.m:` Simple FETI implementation for solving the 2D Helmholtz equation on a square or L-shaped domain using triangular or square finite elements
* `helm_ras.m:` One-level and two-level RAS for solving the 2D Helmholtz equation on a square domain using triangular or square finite elements
* `lshape_ddm2.m:`
* `lshape_ddm.m:`
* `lshape_elasticity.m:`
### Utility Scripts
* `integrals_qfem.m:`
* `integrals_tfem.m:`
### Utility Functions
* `LinearTriangleElementStiffness.m:`
* `lshape_fem_2.m:`
* `lshape_fem.m:`
* `qfem_assemble.m:`
* `qfem.m:`
* `tfem_assemble.m:`
* `tfem.m:`
