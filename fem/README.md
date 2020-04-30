# Finite Element Method (FEM)
### Main Scripts
#### Base FEM Solution Scripts
* `fem3.m:` Solves the **2D Helmholtz** equation on a square domain using **triangular** FEs
* `fem4.m:` Solves the **2D Helmholtz** equation on a square domain using **square** FEs
#### FEM + Preconditioned Krylov Subspace Method Comparisons
* `fem_helmholtz.m:` Benchmarks various **Krylov subspace methods** combined with **Jacobi**, **Gauss-Seidel**, and **Symmetric Gauss-Seidel** preconditioners on the **2D Helmholtz** equation, discretized on a square domain using **triangular** or **square** FEs (see chapter 2 in [doc](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf) for more details)
#### FEM + Domain Decomposition Method on an L-shaped domain
* `lshape_ddm.m:` Solves the **2D Helmholtz** equation on an L-shaped domain using the **Schur complement method**
* `lshape_ddm2.m:` Solves the **2D Helmholtz** equation on an L-shaped domain using the **Schur complement method** combined with **SPAI** and **iLU** preconditioners for the interface problem (see chapter 4 in [doc](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf) for more details)
#### FEM + Advanced Domain Decomposition Methods
* `helm_ras.m:` One-level and two-level **RAS** for solving the **2D Helmholtz** equation on a square domain using **triangular** or **square** FEs
* `helm_feti.m:` Simple **FETI** implementation for solving the **2D Helmholtz** equation on a square or L-shaped domain using **triangular** or **square** FEs
#### Vector FEM for Elasticity + RAS
* `lshape_elasticity.m:` Solves an **elasticity** problem using **vector FEM** and **RAS** with coarse spaces (see chapter 6 in [doc](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf) for more details)
### Utility Scripts
* `tfem_integrals.m:` Symbolically computes the **integrals** for use in FEM with **triangular** elements
* `qfem_integrals.m:` Symbolically computes the **integrals** for use in FEM with **square** elements
### Utility Functions
#### FEM Discretization Functions - Coordinate (COO) and Connectivity (CON) Matrix Computation
* `tfem_discretize.m:` Computes the COO and CON matrices for the discretization of a **square domain** using **triangular** FEs
* `qfem_discretize.m:` Computes the COO and CON matrices for the discretization of a **square domain** using **square** FEs
* `get_lshape_from_square.m:` Computes the COO and CON matrices representing the FEM discretization of an **L-shaped domain** given the discretization of a square domain
* `get_lshape_from_square_2.m:` Computes the COO and CON matrices representing the FEM discretization of an **L-shaped domain** given the discretization of a square domain
#### FEM Linear System Setup Functions
* `tfem_assemble.m:` Sets up the **linear system** given the COO and CON matrices for the case of **triangular** FEs
* `qfem_assemble.m:` Sets up the **linear system** given the COO and CON matrices for the case of **square** FEs
* `LinearTriangleElementStiffness.m:` Computes **stiffness matrices** for solving the **elasticity** problem using **vector FEM**
