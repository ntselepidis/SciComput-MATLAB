# Domain Decomposition Methods / Preconditioners
### Main Scripts
* `ddmeth.m:` Script benchmarking various algebraic one-level and two-level domain decomposition preconditioners
* `schurcomp.m:` Script comparing two preconditioned Schur complement methods
### Utility Functions
#### General DDM Preprocessing Functions
* `inner_outer.m:` Classifies the unknowns of each subdomain into local interior (inner) and local interface (outer) points
* `neighbor_domains.m:` Computes subdomain adjacency matrix / couples of neighboring subdomains
* `vertex_separators.m:` Computes a vertex separator given two adjacent subdomains
* `form_overlap.m:` Computes one-level overlap between adjacent subdomains using BFS
* `form_overlap_ml.m:` Computes multi-level overlap between adjacent subdomains using BFS
* `update_dpnts.m:` Updates inner-outer labeling of subdomain unknowns after computing vertex separators
* `global_separators.m:` Computes a list containing all separators and the associated multiplicities / degrees
* `assemble_B.m:` Computes equality constraint matrix B and its pseudoinverse Bd for use in TI schemes
#### Schur Complement Method Functions
* `sc_prep.m:` Schur complement method preprocessing
* `Smultx.m:` Computes y = S * x without assemblying the Schur complement matrix S
* `Srhs.m:` Computes the right-hand-side vector of the Schur complement linear system 
#### Preconditioner Setup Functions
* `scbj_prep.m:` Block Jacobi preconditioner setup for the Schur complement linear system
* `scbp_prep.m:` Block probing preconditioner setup for the Schur complement linear system
* `scas_prep.m:` Exact and sparsified additive Schwarz preconditioner setup for the Schur complement linear system
#### Preconditioner Application Functions
* `blkjac.m:` Block Jacobi preconditioner
* `mti_prec.m:` Modified Tearing and Interconnecting (MTI) preconditioner
* `rasprec.m:` Restricted additive Schwarz (RAS) preconditioner
* `schuras.m:` Additive Schwarz preconditioner for the Schur complement linear system
#### Coarse Spaces for Two-Level DDMs
* `aggregate.m:` Computes restriction matrix for Nicolaides coarse space
* `shem.m:` Computes restriction matrix for SHEM coarse space
#### Probing Technique
* `probe2.m:` Probing technique
#### Core Utility Functions
* `factorize.m:` Computes sparse LU factorization of all subdomains
* `sp_solve.m:` Performs sparse forward/backward substitution
* `ksm.m:` Wrapper of Krylov subspace methods
