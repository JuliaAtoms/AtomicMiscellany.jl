# AtomicMiscellany

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaatoms.org/AtomicMiscellany.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaatoms.org/AtomicMiscellany.jl/dev)
[![Build Status](https://github.com/JuliaAtoms/AtomicMiscellany.jl/workflows/CI/badge.svg)](https://github.com/JuliaAtoms/AtomicMiscellany.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaAtoms/AtomicMiscellany.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAtoms/AtomicMiscellany.jl)

A gallimaufry of disparate utilities useful in atomic physics that are too small and/or
simple to deserve a package of their own. The package aims to provide a consistent and
generic API for all its functionality, and also to cover the features with comprehensive
suit unit tests.

Package includes:

* **Nuclear models.** Methods for evaluating the densities and potentials for various
  spherically symmetric nuclear models, such as the point nucleus, uniform sphere and uniform shell models, Gaussian nucleus etc.
* **Hydrogenic energies.** Methods for calculating the energies of non-relativistic and
  relativistic (Dirac) hydrogen-like bound states, for both point nuclei and some finite
  nuclear models (currently, only the uniform shell nucleus is supported).

For more details, see [the documentation](https://juliaatoms.org/AtomicMiscellany.jl/dev).
