# sde_structural_identifiability
 
`Julia` and `Mathematica` code used to produce the results in the preprint "Exact identifiability analysis for a class of partially observed near-linear stochastic differential equation models" available on arXiv.

## Folder structure

- The folder `figures` contains Julia code used to reproduce Figures 1, 2, and 3.
- The folder `analysis` contains two subfolders: `linear` and `nonlinear` for structural identifiability analysis of each class of model. 
    - `Julia` code in these directories are used to apply the [`StructuralIdentifiability.jl`](https://docs.sciml.ai/StructuralIdentifiability/stable/) package to produce the structurally identifiable quantities for resultant ODE models.
    - `Mathematica` code is used to perform analytical computations and apply Algorithm 1 to analyse the structural identifiability of (in particular) non-linear SDE models. 