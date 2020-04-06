# Scientific Computing MATLAB Toolbox  
This repository contains a collection of MATLAB scripts implementing various computational methods.  
The folder structure of this repository is the following:
* `ddm:` Domain Decomposition Methods / Preconditioners
* `ddm_time:` Domain Decomposition in Time
* `eigen:` Eigenvalue Computations
* `fdm:` Finite Difference Method
* `fem:` Finite Element Method
* `meshless:` Meshless Method
* `multigrid:` Multigrid Method
* `newton_krylov:` Newton-Krylov Method
* `nystrom:` Nystrom Method
* `pipelined_krylov:` Pipelined Krylov Methods
* `spai:` Sparse Approximate Inverse
### Documentation
In each folder there is a separate README describing the use of every script / function in that folder.  
Also, in the doc [A Study of Advanced Computational Methods.pdf](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf), there are various related analyses and results.  
### Setup
First, install the METIS graph partitioning tool from [here](https://github.com/YingzhouLi/metismex).  
Then, run the `set_path.m` script to add METIS as well as all the folders of this repository in the MATLAB path.  
