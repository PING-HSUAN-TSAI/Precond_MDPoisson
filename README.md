# Precond_MDPoisson
Building up Preconditioning technique such as Overlapping Schwarz and Multi grid in Multi-domain Pseudo spectral Poisson Solver

## Cases
### Iterative solver without preconditioning - CG & GMRES
### Iterative solver with Jacobi smoothing   - CG & GMRES + jacobi smoothing
### Iterative solver with Additive Schwarz  - CG
- Might have to consider Overlapping Schwarz to improve the performance ( to be determined )
### Projection method + Additive Schwarz as a wrapper
### CG & GMRES + Additive Schwarz as a wrapper ( working on )
### Multilevel multigrid ( to be done )
