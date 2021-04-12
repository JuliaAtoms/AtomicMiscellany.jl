```@meta
CurrentModule = AtomicMiscellany
```

# Hydrogenic energies

Hydrogenic bound state energies can be calculated using the [`hydrogenic_energy`](@ref) functions.
Depending on the arguments, the function dispatches on different implementations, calculating the energy for different physical cases.

```@docs
hydrogenic_energy
```

See the [related implementation notes](@ref notes-hydrogenic) for more information and details.

### Examples

```@setup ex
using Plots
using AtomicMiscellany: hydrogenic_energy, NRElectron, DiracElectron, UniformShellNucleus, α
```

As a simple usage example, let's plot the PNC energies for a few ``n`` values, comparing the
relativistic and non-relativistic energies.

```@example ex
plot(
    legend=:bottomleft, size = (800, 400),
    xlabel = "Nuclear charge Z", ylabel = "Energy (mc²)",
    title = "n=1-3 PNC energies",
)
Z = range(1, α^-1, length=501)
plot!(Z, hydrogenic_energy.(NRElectron, Z; n = 1) .* α^2, label="NR", c=1)
plot!(Z, hydrogenic_energy.(NRElectron, Z; n = 2) .* α^2, label=false, c=1)
plot!(Z, hydrogenic_energy.(NRElectron, Z; n = 3) .* α^2, label=false, c=1)
plot!(Z, hydrogenic_energy.(DiracElectron, Z; n = 1, κ = -1) .* α^2, label="Dirac", c=2)
plot!(Z, hydrogenic_energy.(DiracElectron, Z; n = 2, κ = -1) .* α^2, label=false, c=2)
plot!(Z, hydrogenic_energy.(DiracElectron, Z; n = 3, κ = -1) .* α^2, label=false, c=2)
```

Similarly, we can look at the FNC correction to the point nucleus energies

```@example ex
plot(
    legend=:bottomright, yaxis = :log10, size = (800, 400),
    xlabel = "Nuclear charge Z", ylabel = "FNC energy correction (mc²)",
    title = "n=1-4 FNC corrections (relativistic)",
)
Z = range(1, α^-1, length=501)
nm = UniformShellNucleus(1e-4) # RMS: 1e-4 a.u., for all Z
δE(p, Z; kwargs...) = hydrogenic_energy(p, nm, Z; kwargs...) .- hydrogenic_energy(p, Z; kwargs...)
plot!(Z, δE.(DiracElectron, Z; n = 1, κ = -1) .* α^2, label="n=1")
plot!(Z, δE.(DiracElectron, Z; n = 2, κ = -1) .* α^2, label="n=2")
plot!(Z, δE.(DiracElectron, Z; n = 3, κ = -1) .* α^2, label="n=3")
plot!(Z, δE.(DiracElectron, Z; n = 4, κ = -1) .* α^2, label="n=4")
```

## Particles

Different free-particle dynamics (e.g. non-relativistic or relativistic) lead to different energies.
We use subtypes of [`AbstractParticle`](@ref) to dispatch on the different types.

```@docs
AbstractParticle
```

Currently, we support spin-1/2 particles described by either the non-relativistic Schrödinger equation or the relativistic Dirac equation.

```@docs
NRParticle
DiracParticle
```

For convenience, we also instantiate cases with specific masses.

```@docs
NRElectron
DiracElectron
NRMuon
DiracMuon
```

## Index

```@index
Pages = ["hydrogenic.md"]
```
