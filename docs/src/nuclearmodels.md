```@meta
CurrentModule = AtomicMiscellany
```

# Nuclear models

Each nuclear model is implemented as a subtype of [`AbstractNuclearModel`](@ref).
An instance of each nuclear model is assumed to be normalized to ``Z = 1``.

```@docs
AbstractNuclearModel
```

The package implements various nuclear models, illustrated by the following plot.
The root-mean-square radius of the charge density is the same for each mode
to keep the potentials comparable.

```@example
using AtomicMiscellany.NuclearModels, Plots
rs = 10 .^ range(-6, -2, length=501)
plot(
    title = "Nuclear potentials with RMS = 10⁻⁴ a.u.",
    xlabel = "r (a.u.)", ylabel = "Nuclear potential V",
    xaxis=:log10, ylims = (-2.1e4, 0.1e4), legend=:bottomright,
)
V(M) = r -> potential(M, r)
plot!(rs, V(PointNucleus()).(rs), label = "PointNucleus")
plot!(rs, V(from_rms(UniformShellNucleus, 1e-4)).(rs), label = "UniformShellNucleus")
plot!(rs, V(from_rms(UniformSphericalNucleus, 1e-4)).(rs), label = "UniformSphericalNucleus")
plot!(rs, V(from_rms(GaussianNucleus, 1e-4)).(rs), label = "GaussianNucleus")
#plot!(rs, V(from_rms(FermiNucleus, 1e-4)).(rs), label = "FermiNucleus") # hide
plot!() # hide
```

```@docs
PointNucleus
UniformShellNucleus
UniformSphericalNucleus
GaussianNucleus
```

## Root-mean-square radius

The following methods are for determining the RMS radius of nuclei, loosely based on
experimental data.

```@docs
rms(::Type{JohnsonSoff1985}, ::Real)
```

## Utilities

```@docs
potential
density
rms
from_rms
AtomicMiscellany.NuclearModels
```

## References

```@docs
Andrae2000
JohnsonSoff1985
```

## Index

```@index
Pages = ["nuclearmodels.md"]
```
