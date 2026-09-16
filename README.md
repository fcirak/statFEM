# statFEM

statFEM is CSMLab’s in-house finite element library, developed over several years by many brilliant [PhD students and postdocs](https://www.csmlab.org/people/) to support research in computational mechanics. The C++ library provides the computational framework underlying a range of our work on finite element methods, solid and structural mechanics, subdivision and isogeometric methods, optimisation, and probabilistic computational mechanics.

We are releasing the library publicly to make this research software available to the wider computational mechanics community.

## Getting started

- **Video series**: [Beginner's guide to statFEM and openFTL](https://www.youtube.com/playlist?list=PLSMyGJmbhFw0QlQvSiQgnUCeW7R7gNaH6) — walks through installing and running the code.
- **Wiki**: [github.com/fcirak/statFEM/wiki](https://github.com/fcirak/statFEM/wiki) — install steps, VS Code setup, Eigen solver notes (Intel Pardiso / Apple Accelerate), and contribution conventions.
- See [Building and running](#building-and-running) below for a quick-start; the wiki has the full, platform-specific version.

## Publications

This code has been used to produce results for:

- Girolami, M., Febrianto, E., Yin, G., Cirak, F. (2021). [The statistical finite element method (statFEM) for coherent synthesis of observation data and model predictions](https://doi.org/10.1016/j.cma.2020.113533). *Computer Methods in Applied Mechanics and Engineering*, 375, 113533.
- Febrianto, E., Butler, L., Girolami, M., Cirak, F. (2022). [Digital twinning of self-sensing structures using the statistical finite element method](https://doi.org/10.1017/dce.2022.28). *Data-Centric Engineering*, 5, e29.
- Koh, K. J., Cirak, F. (2023). [Stochastic PDE representation of random fields for large-scale Gaussian process regression and statistical finite element analysis](https://www.sciencedirect.com/science/article/abs/pii/S0045782523004826). *Computer Methods in Applied Mechanics and Engineering*, 417 Part B, 116353. Poisson and thin-shell examples in this paper correspond to the `del2` and `gshell` modules respectively.

If you use this code, please cite the 2021 paper above.

## Repository structure

| Path | Contents |
|---|---|
| `corlib/` | Core finite element library: meshes, elements, assemblers, linear solvers. Shared by every physics module below. |
| `statX/` | The statistical layer: Gaussian process regression, MCMC, and the statFEM data/model-synthesis machinery itself. |
| `del2/` | Laplace/Poisson problems (∇²). Simplest module — a good place to start. |
| `beam/` | 1D beam elements. |
| `solid/` | 3D solid mechanics (static and dynamic elasticity). |
| `gshell/` | Thin shell elements, built on subdivision surfaces. |
| `subdiv/` | Subdivision-surface geometry and meshing utilities used by `gshell`. |
| `tools/` | Mesh format converters, mesh partitioners, and vendored third-party libraries (NLopt, BDDCML, etc.). |
| `config/` | Per-machine build configuration (see below). |
| `docu/` | Design documentation. |

Each physics module follows the same pattern:
- `apps/reference/` — curated example applications, several of which reproduce results from the
  papers above. Treat these as read-only; copy one into `apps/user/` before modifying it.
- `apps/user/` — individual researchers' applications, at varying stages of completion and
  maintenance. Some are experimental or unfinished.

## Building and running

### Dependencies

- A C++17 compiler (`g++`) and an MPI C++ wrapper (`mpicxx`)
- [Boost](https://www.boost.org/)
- [Eigen](https://eigen.tuxfamily.org/)
- [PETSc](https://petsc.org/) (for the parallel/PETSc-based solid solvers)
- [SuperLU](https://portal.nersc.gov/project/sparse/superlu/)
- [GSL](https://www.gnu.org/software/gsl/)
- [Scotch](https://www.labri.fr/perso/pelegrin/scotch/) and [METIS](https://github.com/KarypisLab/METIS)/ParMETIS (mesh partitioning)
- LAPACK/BLAS, and a Fortran runtime (`gfortran`) for linking Fortran-based solver libraries
- [gmsh](https://gmsh.info/), for generating meshes some example applications expect
- Python 3.9–3.12 via [Poetry](https://python-poetry.org/), for the `baryrat`-based rational
  approximation used by `statX`'s rational SPDE precision-matrix code (see note below)

### Environment and machine configuration

1. Set the `STATFEMROOT` environment variable to the path of this repository, e.g. add
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

Each application directory has its own `Makefile`; there is no repo-wide build step beyond
setting `STATFEMROOT` and picking a per-app directory to build in.

### VS Code

The repository also works with VS Code's Makefile Tools extension. The root `Makefile`'s
`SUBDIRS` variable points at a single app directory at a time — edit it by hand to switch which
app the extension's build/launch targets operate on. See the wiki's
[development and usage guide](https://github.com/fcirak/statFEM/wiki/statFEM-development-and-usage)
for the full setup, including IntelliSense configuration.

## Contributing

See the wiki's [GitHub best-practice guide](https://github.com/fcirak/statFEM/wiki/Github-statFEM-best-practice)
for branch naming, issue, and pull-request conventions used in this repository.

## License

Not yet specified.
