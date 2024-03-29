# Scientific Computing MATLAB Toolbox  
This repository contains a collection of MATLAB scripts implementing various computational methods.  
The folder structure of this repository is the following:
```
SciComput-MATLAB
├── ddm                  <-- Domain Decomposition Methods / Preconditioners
├── ddm_time             <-- Domain Decomposition in Time
├── eigen                <-- Eigenvalue Computations
├── fdm                  <-- Finite Difference Method
├── fem                  <-- Finite Element Method
├── integration          <-- Numerical Integration Methods
├── krylov               <-- Krylov Subspace Methods
├── low_rank_approx      <-- Low-Rank Matrix Approximation Methods
├── meshless             <-- Meshless Method
├── multigrid            <-- Multigrid Method
├── newton_krylov        <-- Newton-Krylov Method
├── nystrom              <-- Nystrom Method
└── spai                 <-- Sparse Approximate Inverse
```
### Documentation
* In each folder there is a separate `README` describing the use of every script / function in that folder.  
* Also, in the doc [A Study of Advanced Computational Methods.pdf](https://github.com/ntselepidis/SciComput-MATLAB/blob/master/A%20Study%20of%20Advanced%20Computational%20Methods.pdf), there are various related analyses and results.  
### Installation
* First, install the METIS graph partitioning tool from [here](https://github.com/YingzhouLi/metismex).  
* Then, download the Scientific Computing Toolbox, by running the command:
```bash
git clone https://github.com/ntselepidis/SciComput-MATLAB.git # Download Scientific Computing Toolbox
``` 
* And finally, set the MATLAB path by executing the `set_path.m` script.
