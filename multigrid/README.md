# Multigrid Method
### Main Scripts
#### Algebraic Multigrid (AMG)
* `amg2d.m:` V/W-cycle algebraic multigrid with **damped Jacobi** smoother solving **2D Poisson** equation
* `amg2d_spai.m:` V/W-cycle algebraic multigrid with **SPAI** smoother solving **2D Poisson** equation
#### Geometric Multigrid (GMG)
* `multigrid2d_convection_diffusion.m:` V/W-cycle geometric multigrid with **damped Jacobi/Gauss-Seidel** smoother solving **2D convection-diffusion** equation (see chapter 1 in [doc](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf) for more details)
* `multigrid2d.m:` V/W-cycle geometric multigrid with **damped Jacobi** smoother solving **2D Poisson** equation
* `multigrid2d_spai.m:` V/W-cycle geometric multigrid with **SPAI** smoother solving **2D Poisson** equation
* `multigrid.m:` V/W-cycle geometric multigrid with **damped Jacobi** smoother solving **1D Poisson** equation
### Utility Functions
* `coarsening.m:` Coarsening algorithm for algebraic multigrid
* `mgmu2.m:` Multigrid cycle with **explicit** smoothing
* `mgmu3.m:` Multigrid cycle with **implicit** smoothing
* `mgmu.m:` Multigrid cycle with **damped Jacobi** smoothing

