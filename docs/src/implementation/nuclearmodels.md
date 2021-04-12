```@meta
CurrentModule = AtomicMiscellany
```

# [Notes: nuclear models](@id notes-nuclei)

Notes on the internal implementation details of the nuclear models. Anything discussed here
is **not** considered to be part of the public API.

## Density vs. potential

The nuclear charge density and the corresponding potential are related to each other via the
Poisson equation. In radial coordinates, for a spherically symmetric density, the potential
at ``r``, as per equation (22) in [`Andrae2000`](@ref), is given by the following integral

```math
V(r) = - \frac{4\pi}{r}\left[
    \int_0^r s^2 \rho(s) ~ds
    + r \int_r^\infty s \rho(s) ~ds
\right]
```

In the tests, we use this to numerically verify that the implementations for [`density`](@ref)
and [`potential`](@ref) actually match, as they are generally implemented independently
based on analytical expressions.

The numerical integration can be done using the following routine

```@example integrate
using QuadGK
function potential_quadk(ρ, r; rmax = Inf, rscale = 1.0)
    # Based on equation (22) in Andrae. The potential is split into two separate integrals
    # on the intervals [0, r] and [r, ∞).
    #
    # However, often the density can have its non-trivial parts at very small r values,
    # which means that the adaptive quadgk algorithm might not realize that it even has
    # non-zero parts. To work around that, we do 2 things:
    #
    # 1. We scale the r-axis with rscale (so r = rscale*t), and actually integrate in t.
    #    In practice, we'll use the RMS value for the scaling, since we'd expect that to
    #    correspond to the scale at which the density has a non-trivial shape.
    # 2. If r is much greater than rscale, we split the [0, r] integration into two
    #    sub-intervals, hoping that we capture all the non-trivial part in the first
    #    interval.
    #
    # The change of variable r = rscale*t leads to
    #
    #   ∫ f(r) dr = rscale ∫ f(rscale*t) dt
    #
    # So we have to also multiply the resulting integral with rscale.
    _I1(t) = rscale^3 * t^2 * ρ(rscale*t)
    I1, δ1 = min(r, rmax)/rscale > 10.0 ? quadgk(_I1, 0, 10.0, min(r, rmax)/rscale) :
        quadgk(_I1, 0, min(r, rmax)/rscale)
    I2, δ2 = (r < rmax) ? quadgk(t -> rscale^2 * t * ρ(rscale*t), r/rscale, rmax/rscale) : (0.0, 0.0)
    # add up the subintegral based on equation (22):
    -4π*(I1/r + I2)
end
nothing # hide
```

As an example, let's plot a few arbitrary cases

```@example integrate
using AtomicMiscellany.NuclearModels
using Plots: plot, plot!
rs = range(0, 1, length=101)
plot(legend=:bottomright)
let nm = UniformSphericalNucleus(0.6)
    plot!(
        rs, potential_quadk.(r -> density(nm, r), rs, rmax = 0.6, rscale = rms(nm)),
        label = "Uniform sphere", c = 1, lw = 3,
    )
    plot!(rs, (r -> potential(nm, r)).(rs), c=:black, ls=:dash, label=nothing)
end
let nm = GaussianNucleus(0.25)
    plot!(
        rs, potential_quadk.(r -> density(nm, r), rs, rscale = rms(nm)),
        label = "Gaussian", c = 2, lw = 3,
    )
    plot!(rs, (r -> potential(nm, r)).(rs), c=:black, ls=:dash, label="Analytic")
end
```
