# Multigrid Method
### Main Scripts
#### Algebraic Multigrid (AMG)
* `amg2d.m:` V/W-cycle AMG with **damped Jacobi** smoother solving **2D Poisson** equation
* `amg2d_spai.m:` V/W-cycle AMG with **SPAI** smoother solving **2D Poisson** equation
#### Geometric Multigrid (GMG)
* `multigrid.m:` V/W-cycle GMG with **damped Jacobi** smoother solving **1D Poisson** equation
* `multigrid2d.m:` V/W-cycle GMG with **damped Jacobi** smoother solving **2D Poisson** equation
* `multigrid2d_spai.m:` V/W-cycle GMG with **SPAI** smoother solving **2D Poisson** equation
* `multigrid2d_convection_diffusion.m:` V/W-cycle GMG with **damped Jacobi/Gauss-Seidel** smoother solving **2D convection-diffusion** equation (see chapter 1 in [doc](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf) for more details)
### Utility Functions
* `mgmu.m:` Core mu-cycle MG with **damped Jacobi** smoothing
* `mgmu2.m:` Core mu-cycle MG with **explicit** smoothing
* `mgmu3.m:` Core mu-cycle MG with **implicit** smoothing
* `coarsening.m:` Coarsening algorithm for AMG

