# Finite Element Method (FEM)
### Main Scripts
#### Base FEM Solution Scripts
* `fem3.m:` Solves **2D Helmholtz** equation on a square domain using **triangular** finite elements
* `fem4.m:` Solves **2D Helmholtz** equation on a square domain using **square** finite elements
#### FEM + Preconditioned Krylov Subspace Method Comparisons
* `fem_helmholtz.m:` Benchmarks various **Krylov subspace methods** combined with **Jacobi**, **Gauss-Seidel**, and **Symmetric Gauss-Seidel** preconditioners on solving the **2D Helmholtz** equation on a square domain using **triangular** or **square** finite elements (see chapter 2 in [doc](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf) for more details)
#### FEM + Domain Decomposition Method on an L-shaped domain
* `lshape_ddm.m:` Solves **2D Helmholtz** equation on an L-shaped domain using the **Schur complement method**
* `lshape_ddm2.m:` Solves **2D Helmholtz** equation on an L-shaped domain using the **Schur complement method** combined with **SPAI** and **ILU** preconditioners for the interface problem (see chapter 4 in [doc](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf) for more details)
#### FEM + Advanced Domain Decomposition Methods
* `helm_ras.m:` One-level and two-level **RAS** for solving the **2D Helmholtz** equation on a square domain using **triangular** or **square** finite elements
* `helm_feti.m:` Simple **FETI** implementation for solving the **2D Helmholtz** equation on a square or L-shaped domain using **triangular** or **square** finite elements
#### Vector FEM for Elasticity + RAS
* `lshape_elasticity.m:` Solves **elasticity** problem using **vector FEM** and **RAS** with coarse spaces (see chapter 6 in [doc](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf) for more details)
### Utility Scripts
* `integrals_tfem.m:` Symbolically computes the **integrals** for use in FEM with **triangular** elements
* `integrals_qfem.m:` Symbolically computes the **integrals** for use in FEM with **square** elements
### Utility Functions
* `tfem.m:` Computes the **coordinate** and **connectivity** matrix for the discretization of a **square domain** using **triangular** finite elements
* `qfem.m:` Computes the **coordinate** and **connectivity** matrix for the discretization of a **square domain** using **square** finite elements
* `lshape_fem.m:` Computes the **coordinate** and **connectivity** matrix representing the FEM discretization of an **L-shaped domain** given the discretization of a square domain
* `lshape_fem_2.m:` Computes the **coordinate** and **connectivity** matrix representing the FEM discretization of an **L-shaped domain** given the discretization of a square domain
* `tfem_assemble.m:` Assembles the **linear system** given the coordinate and connectivity matrices for the case of **triangular** finite elements
* `qfem_assemble.m:` Assembles the **linear system** given the coordinate and connectivity matrices for the case of **square** finite elements
* `LinearTriangleElementStiffness.m:` Computes **stiffness matrices** for solving the **elasticity** problem using **vector FEM**
