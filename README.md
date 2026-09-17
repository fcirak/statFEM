# statFEM

statFEM provides an implementation of the statistical finite element method. It is based on [CSMLab](https://www.csmlab.org/)’s in-house finite element library, developed over several years by many brilliant [PhD students and postdocs](https://www.csmlab.org/people.html). The original C++ library forms the computational framework underlying a range of our work on finite element methods, solid and structural mechanics, subdivision and isogeometric methods, optimisation, and probabilistic computational mechanics.

![Digital Twin](https://raw.githubusercontent.com/wiki/fcirak/statFEM/images/intro_snapshot.png)

## Relevant publications

- Girolami, M., Febrianto, E., Yin, G., Cirak, F. (2021). [The statistical finite element method (statFEM) for coherent synthesis of observation data and model predictions](https://doi.org/10.1016/j.cma.2020.113533). *Computer Methods in Applied Mechanics and Engineering*, 375, 113533.
- Febrianto, E., Butler, L., Girolami, M., Cirak, F. (2022). [Digital twinning of self-sensing structures using the statistical finite element method](https://doi.org/10.1017/dce.2022.28). *Data-Centric Engineering*, 3, e31.
- Koh, K. J., Cirak, F. (2023). [Stochastic PDE representation of random fields for large-scale Gaussian process regression and statistical finite element analysis](https://doi.org/10.1016/j.cma.2023.116358). *Computer Methods in Applied Mechanics and Engineering*, 417, Part B, 116358. 

See [Google Scholar](https://scholar.google.com/citations?hl=en&user=01BvTAsAAAAJ&view_op=list_works&sortby=pubdate) for further publications, with the most recent papers first.

## Repository structure

| Path | Contents |
| --- | --- |
| `corlib/` | Core finite element library: meshes, elements, assemblers, linear solvers. Shared by every physics module below. |
| `statX/` | Gaussian process regression, MCMC, statFEM helpers. |
| `del2/` | Laplace/Poisson problems. |
| `beam/` | Beam problems. |
| `solid/` | Solid problems. |
| `gshell/` | Isogeometric subdivision shell problems. |
| `subdiv/` | Subdivision surfaces and basis functions. |
| `tools/` | Mesh format converters, mesh partitioners, and third-party libraries (NLopt, etc.). |
| `config/` | Machine-specific configuration files (see below). |
| `docu/` | Documentation. |

## Building and running

### Getting started

- **Video series**: [Beginner's guide to statFEM and openFTL](https://www.youtube.com/playlist?list=PLSMyGJmbhFw0QlQvSiQgnUCeW7R7gNaH6) — walks through installing and running the code.
- **Wiki**: [github.com/fcirak/statFEM/wiki](https://github.com/fcirak/statFEM/wiki) — install steps, VS Code setup, Eigen solver notes (Intel Pardiso / Apple Accelerate), and contribution conventions.

### Dependencies

- [Boost](https://www.boost.org/)
- [Eigen](https://eigen.tuxfamily.org/)
- [PETSc](https://petsc.org/) (for the parallel/PETSc-based solid solvers)
- [SuperLU](https://portal.nersc.gov/project/sparse/superlu/)
- [METIS](https://github.com/KarypisLab/METIS)/ParMETIS (mesh partitioning)
- LAPACK/BLAS and a Fortran runtime (`gfortran`) for linking Fortran-based solver libraries
- [Gmsh](https://gmsh.info/) for generating meshes required by some example applications

### Environment and machine configuration

1. Set the `STATFEMROOT` environment variable to the path of this repository, e.g. by adding
   `export STATFEMROOT=/path/to/statFEM` to your shell profile.
2. The build reads compiler/library paths from `config/statFEM_definitions.<hostname>`, selected
   automatically by your machine's hostname (`hostname -s`), falling back to
   `config/statFEM_definitions.default` if no host-specific file exists. If the defaults don't
   match your system, copy `config/statFEM_definitions.default` to
   `config/statFEM_definitions.<your-hostname>` and adjust the paths.

### Building and running an application

```sh
cd del2/apps/reference/laplace
make
./laplace
```

Each application directory has its own `Makefile`. There is no repo-wide build step beyond
setting `STATFEMROOT`.

### VS Code

The repository also works with VS Code's Makefile Tools extension. The root `Makefile`'s
`SUBDIRS` variable points at a single app directory at a time — edit it by hand to switch which
app the extension's build/launch targets operate on. See the wiki's
[development and usage guide](https://github.com/fcirak/statFEM/wiki/statFEM-development-and-usage)
for the full setup, including IntelliSense configuration.

## Getting in touch

If you have any questions, feedback, or suggestions, please feel free to contact us at [statfem@csmlab.org](mailto:statfem@csmlab.org).
