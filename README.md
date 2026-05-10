# GLOM Hamiltonian Structure: Code Repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20112613.svg)](https://doi.org/10.5281/zenodo.20112613)

**Paper:** Quadratic invariants and Hamiltonian structure in coupled gyrostat low-order model hierarchies  
**Authors:** Seshadri and Lakshmivarahan (2026a)  
**Journal:** (submitted)

This repository contains MATLAB code supporting the computations in the paper and its
Supplementary Information (SI). All scripts require the **MATLAB Symbolic Math Toolbox**.

---

## Repository structure

```
GLOM_Hamiltonian_repo/
│
├── invariants/          # Quadratic invariant computation
│   ├── colechelon_symbolic_single_gyrostat.m
│   ├── colechelon_symbolic_model2.m
│   ├── colechelon_symbolic_model3.m
│   ├── getrank_model1_all_subclasses.m
│   └── getrank_model2_all_subclasses.m
│
├── hamiltonian/         # Jacobi condition for canonical hierarchy
│   ├── Jacobi_sparse_hierarchy.m
│   └── Jacobi_dense_hierarchy.m
│
├── hierarchy/           # Casimir functions and gradient consistency
│   ├── Casimir_sparse_hierarchy.m
│   ├── Casimir_hierarchy_of_model4.m
│   ├── Casimir_hierarchy_of_model5.m
│   └── Casimir_gradient_consistency.m
│
└── 3K_exhaustive/       # Exhaustive 3^K representation search
    ├── Jacobi_exhaustive_model1.m
    ├── Jacobi_exhaustive_model2.m
    └── Jacobi_exhaustive_model2_K3.m
```

---

## Folder descriptions

### `invariants/` — Standard algebraic approach

These scripts compute quadratic invariants using the constraint-matrix (column-echelon)
method of the SI. For a GLOM with $M$ modes, the invariant $C = \sum_i d_i x_i^2/2 + \sum_{i < j} e_{ij} x_i x_j + \sum_i f_i x_i$ must satisfy $\dot{C} = 0$, which gives a linear system $Au = 0$ whose null-space dimension
equals the invariant count.

| File | What it computes |
|---|---|
| `colechelon_symbolic_single_gyrostat.m` | Invariants of a single gyrostat ($K=1$, $M=3$) |
| `colechelon_symbolic_model2.m` | Invariants of Model 2 ($K=2$, $M=5$, sparse) |
| `colechelon_symbolic_model3.m` | Invariants of Model 3 ($K=3$, $M=5$, dense) |
| `getrank_model1_all_subclasses.m` | All $2^{12}$ subclasses of Model 1: invariant count, Hamiltonian flag, rank($J$), Casimir |
| `getrank_model2_all_subclasses.m` | All $2^{12}$ subclasses of Model 2 (same outputs as above) |

### `hamiltonian/` — Jacobi condition for the canonical $L_A$ representation

These scripts build $J = \sum_k J^{(k)}$ using the canonical (type-$A$) representation
for each gyrostat and evaluate the Jacobi identity symbolically at each hierarchy level.

| File | What it computes |
|---|---|
| `Jacobi_sparse_hierarchy.m` | Jacobi condition at $K=1,2,3,4$ for sparse hierarchy; incremental condition per level |
| `Jacobi_dense_hierarchy.m` | Jacobi condition at $K=1,\ldots,6$ for dense hierarchy |

### `hierarchy/` — Casimir functions and gradient consistency

These scripts compute Casimir functions explicitly and verify the gradient consistency
property: the $K$-level Casimir gradient is collinear with the projection of the
$(K+1)$-level gradient.

| File | What it computes |
|---|---|
| `Casimir_sparse_hierarchy.m` | Casimirs at $K=1,2,3,4$ for sparse hierarchy with $q_k=0$ |
| `Casimir_hierarchy_of_model4.m` | Casimirs for 2D Rayleigh–Bénard hierarchy (Model 4) |
| `Casimir_hierarchy_of_model5.m` | Casimirs for 3D Rayleigh–Bénard hierarchy (Model 5) |
| `Casimir_gradient_consistency.m` | Null-space intersection computation verifying Theorem 2 (Appendix 4 of main text) |

### `3K_exhaustive/` — Exhaustive $3^K$ representation search

Each gyrostat's contribution to $J$ can be written in three ways ($L_A$, $L_B$, $L_C$),
giving $3^K$ candidate Poisson matrices for $K$ gyrostats. These scripts compute the
Jacobi condition for every candidate symbolically.

| File | What it computes |
|---|---|
| `Jacobi_exhaustive_model1.m` | All $3^2=9$ representations of Model 1 |
| `Jacobi_exhaustive_model2.m` | All $3^2=9$ representations of Model 2; identifies the unconditionally Hamiltonian $(L_A, L_C)$ case |
| `Jacobi_exhaustive_model2_K3.m` | All $3^3=27$ representations of the $K=3$ sparse hierarchy |

---

## Requirements

- MATLAB R2019b or later (earlier versions are untested)
- Symbolic Math Toolbox (required for all scripts)

---

## Citation

If you use this code, please cite:

> Seshadri, A.K. and Lakshmivarahan, S. (2026a). Quadratic invariants and Hamiltonian
> structure in coupled gyrostat low-order model hierarchies. [doi pending]
>
> Code: akseshadri. (2026). akseshadri/GLOM-Hamiltonian: Initial release (v1.0.0).
> Zenodo. https://doi.org/10.5281/zenodo.20112613

---

## Correspondence

For questions about the code or paper, contact **ashwins@iisc.ac.in** or open a GitHub issue.
