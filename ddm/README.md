# Domain Decomposition Methods / Preconditioners
### Main Scripts
* `ddmeth.m:`
* `schurcomp.m:`
### Utility Functions
#### General Domain Decomposition Preprocessing Functions
* `inner_outer.m:`
* `neighbor_domains.m:`
* `vertex_separators.m:`
* `form_overlap.m:`
* `form_overlap_ml.m:`
* `update_dpnts.m:`
* `global_separators.m:`
* `assemble_B.m:`
#### Schur Complement Method Functions
* `sc_prep.m:`
* `Smultx.m:`
* `Srhs.m:`
#### Preconditioner Setup Functions
* `scbj_prep.m:`
* `scbp_prep.m:`
* `scas_prep.m:`
#### Preconditioner Application Functions
* `blkjac.m:`
* `mti_prec.m:`
* `rasprec.m:`
* `schuras.m:`
#### Coarse Spaces
* `aggregate.m:`
* `shem.m:`
#### Probing Technique
* `probe2.m:`
#### Core Utility Functions
* `factorize.m:`
* `sp_solve.m:`
* `ksm.m:`
