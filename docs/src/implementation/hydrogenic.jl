# ```@meta
# CurrentModule = AtomicMiscellany
# ```
#
# # [Notes: hydrogenic energies](@id notes-hydrogenic)

using Pkg #src
Pkg.activate(joinpath(@__DIR__, "..", "..")) #src
using AtomicMiscellany: AtomicMiscellany, hydrogenic_energy, UniformShellNucleus, NRElectron, DiracElectron, α #hide
using Plots, Formatting #hide

# ## Point nucleus
#
# The energies for point nuclei (represented with [`PointNucleus`](@ref)) are calculated using analytic formulae.
#
# For a non-relativisic particles ([`NRParticle`](@ref)) the energy is given by
#
# ```math
# E_n = - \frac{(Z\alpha)^2}{2 n^2} mc^2
# ```
#
# and for relativistic particles ([`DiracParticle`](@ref)) by
#
# ```math
# m c^2 \left\{\left[1 + \left(
# \frac{Z\alpha}{n - |\kappa| + \sqrt{\kappa^2 - (Z\alpha)^2}}
# \right)^2 \right]^{-1/2} - 1 \right\}
# ```
#
# The only thing to note is that for the relativistic case, we set the ``E = 0`` for ``Z = 0``, so that the convention for the energy scale for the relativistic and non-relativistic cases would be the same.

# ## Uniform spherical shell
#
# There are no analytical expressions for the energies for finite nuclei. However, in the case of the uniform spherical shell ([`UniformShellNucleus`](@ref); aka. cut-off Coulomb, top-slice model), the potential is simple and the equations solvable in the ranges ``[0, R]`` and ``[R, \infty)`` separately.
#
# Following the expressions in [`Andrae2000`](@ref), we implement a root-finding based algorithm for determining the energies, in both the relativistic and non-relativistic case.

## We use this function below to showcase the performance of the implementation of hydrogenic_energies.
function plot_fnc_energies(p; n,
    Rs = [0.5, 1e-1, 5e-2, 2e-2, 1e-2, 7e-3, 3e-3, 1e-3, 1e-4],
    qnumbers...
)
    Z = range(1, isa(p, AtomicMiscellany.DiracParticle) ? α^-1 : 150, length=1001)
    emin = hydrogenic_energy.(p, Z; n = n, qnumbers...)
    emax = hydrogenic_energy.(p, Z; n = n + 1, qnumbers...)
    p1 = plot(hcat(Z,Z), hcat(emin,emax), c=:black, ls=:dash, legend=:bottomleft, label=false)
    p2 = plot(hcat(Z,Z), hcat(emin ./ abs.(emin), emax ./ abs.(emin)), label=false, c=:black, ls=:dash)
    for (i, R) = enumerate(Rs)
        es = map(Z) do Z
            try
                hydrogenic_energy(p, UniformShellNucleus(R), Z; n = n, qnumbers...)
            catch e
                return NaN
            end
        end
        plot!(p1, Z, es, label="R=$R", c = i)
        plot!(p2, Z, es./ abs.(emin), label=false, c = i)
    end
    ylabel!(p1, "E"); ylabel!(p2, "E / E(PNC)")
    xlabel!(p2, "Nuclear charge Z")
    plot(p1, p2, layout=(2,1), size=(800, 600))
    vline!([α^-1 α^-1], label=false, c=:black, ls=:dot)
    title!("n = $n, $(qnumbers...)", subplot=1)
end
nothing # hide

# ### Non-relativisic energies
#
# Calculates the matching function ``f`` derived from [`Andrae2000`](@ref), equations (247-248, 256-258). ``f`` is derived by multiplying the equation through with the denominators to remove any infinities these denominators might bring.
#
# The function `AtomicMiscellany.f_andrae_nonrel_topslice` implements the following function, which is then fed into a root-finder:
#
# ```math
# \begin{aligned}
# f(E, r, Z, R, \ell) &=
# \beta j_\ell(x_{\rm{in}}) U(a, b, x_{\rm{out}})
# \\ &+ 2 \beta a j_\ell(x_{\rm{in}}) U(a+1, b+1, x_{\rm{out}})
# \\ &- \alpha j_{\ell+1}(x_{\rm{in}}) U(a, b, x_{\rm{out}}) ≡ 0
# \end{aligned}
# \\[1em]
# \alpha_{n\ell} = \sqrt{2(E + Z/R)}, \quad
# x_{\rm{in}} = \alpha_{n\ell} r, \quad
# \beta_{n\ell} = \sqrt{-2E}, \quad
# x_{\rm{out}} = 2 \beta_{n\ell} r,
# \\[1em]
# a = \frac{b}{2} - \nu_{n\ell}, \quad
# b = 2(\ell + 1), \quad
# \nu_{n\ell} = Z/\beta_{n\ell}
# ```

## We can use this function to visualize f:
function plot_nonrel_fnc_f(; R, Z, ℓ)
    emin = hydrogenic_energy(NRElectron, Z; n = 1)
    ## Excluding the first point, since at E = -Z/R, f is zero, which causes problems
    es = range(max(1.1emin, -Z/R), 0, length=502)[2:end-1]
    ys = (E -> AtomicMiscellany.f_andrae_nonrel_topslice(Z=Z, r=R, R=R, ℓ=ℓ, E = E)).(es)
    any(iszero.(ys)) && @warn "Zeros in ys"
    plot(es, abs.(ys), yaxis=:log10, legend=false, c = (y -> y>0 ? 2 : 1).(ys))
    absy_sorted = filter(!iszero, sort(abs.(ys)))
    ylims!(absy_sorted[1]/2, absy_sorted[490])
    e0s = map(ℓ .+ (1:6)) do n
        try
            hydrogenic_energy(NRElectron, UniformShellNucleus(R), Z; n = n, ℓ = ℓ, verbose=false)
        catch
            NaN
        end
    end
    xlabel!("E"); ylabel!("log |f|")
    title!("Z=$(fmt(".2f", Z)), ℓ=$ℓ, R=$(fmt(".2e", R))")
    vline!(e0s, c=:gray)
    vline!([hydrogenic_energy(NRElectron, Z; n = n) for n = 1:3], c=:black, ls=:dash)
    plot!(size=(800, 400))
end
nothing # hide

# Here is ``f(E)`` for ``Z = 1`` and ``\ell = 0``, for two very different nuclear sizes.
# ``R = 10^{-4}~\mathrm{a.u.}`` is approximately physical, whereas ``R = 2.5~\mathrm{a.u.}`` is unphysically huge. However, the latter is useful for showing how the graph of the function shifts as we make the nucleus bigger and bigger.
#
# The black dashed vertical lines indicate the hydrogenic energies for point nuclei, which can be determined using explicit formulae. As we need to specify a bracketing interval for the root-finder, we use those energies.
# This seems like a reasonable choice, as we know that the FNC energy will always we higher than the corresponding PNC energy, but unlikely to reach the PNC energy of the 2s, unless the parameters are very extreme.

plot(
    plot_nonrel_fnc_f(R = 1e-4, Z = 1, ℓ = 0),
    plot_nonrel_fnc_f(R = 2.5, Z = 1, ℓ = 0),
    layout = (2, 1), size=(800, 600),
)

# To illustrate the potential problems with the bracketing intervals, let's see what happens (for a huge nucleus) when we increase ``Z`` a lot. Eventually, we will hit the `2s` energy. However, this should not really happen in any physically relevant cases.

a = @animate for Z = range(1, 550, length=101)
    plot_nonrel_fnc_f(R = 1e-2, Z = Z, ℓ = 0)
end
gif(a, fps=10)

# To showcase the algorithm an its limits, here are the `1s` energies for various nuclear sizes (`R` values).
# The black dashed lines once again indicate the bracketing interval of PNC energies.

plot_fnc_energies(NRElectron, n = 1, ℓ = 0)

# Now, if we do the same thing for ``n = 2``, we start seeing weird cases where

plot_fnc_energies(NRElectron, n = 2, ℓ = 0)

# This can easily be understood by looking at the behaviour of ``f`` as you increase ``Z`` for a large nucleus (`R = 0.1`):

a = @animate for Z = range(25, 100, length=101)
    plot_nonrel_fnc_f(R = 0.1, Z = Z, ℓ = 0)
end
gif(a, fps=10)

# At one point, the zero corresponding to the `1s` energy enters the bracketing interval of the `2s`. Root-finder, of course, can not distiguish that, so it just returns the `1s` energy. Again though, this is only an issue if the value of ``R`` or ``Z`` gets extremely large.

# We see the same happening for higher ``\ell`` values too. One thing to note is that we also get in trouble for the lowest eigenvalue, because the ``E = -Z/R``, below which ``f(E)`` becomes complex, enters the ``n=2`` bracketing interval.

plot_fnc_energies(NRElectron, n = 2, ℓ = 1)

#-

a = @animate for Z = range(25, 100, length=101)
    plot_nonrel_fnc_f(R = 0.1, Z = Z, ℓ = 1)
end
gif(a, fps=10)

# However, again, for reasonable ``R`` and ``Z`` values, everything behaves reasonably.
# Also, correctly, we do not get any zeros at lower ``n`` values for higher ``\ell`` values:

plot(
    plot_nonrel_fnc_f(R = 1e-4, Z = 1, ℓ = 0),
    plot_nonrel_fnc_f(R = 1e-4, Z = 1, ℓ = 1),
    plot_nonrel_fnc_f(R = 1e-4, Z = 1, ℓ = 2),
    layout = (3, 1), size=(800, 900),
)

# ### Relativistic energies
#
# The relativistic energy is determined in essentially the same way as in the non-relativistic case.
# The relevant equation in [`Andrae2000`](@ref) is (284), setting up the matching condition between the inner and outer parts of the orbital:
#
# ```math
# L^{\mathrm{in}}_{n\kappa}(R) - \frac{\gamma}{R} = - T^{\mathrm{out}}_{n\kappa}(R)
# ```
#
# ``L^{\mathrm{in}}_{n\kappa}(r)`` is the logarithmic derivative of the ``P`` component of the inner part of the wavefunction (i.e. for ``r \leq R``).
# According to equations (273) and (274), it can be written as
#
# ```math
# L^{\mathrm{in}}_{n\kappa}(r) = \frac{\ell + 1}{r} - T^{\mathrm{in}}_{n\kappa}(r)
# = \frac{\ell + 1}{r} - \alpha_{n\kappa} \frac{j_{\ell+1}(x_{\mathrm{in}})}{j_{\ell}(x_{\mathrm{in}})}
# \\[1em]
# \alpha_{n\kappa} = \sqrt{2 (E + Z/R) + [(E + Z/R)/c]^2},
# \quad
# x_{\mathrm{in}} = \alpha_{n\kappa} r
# ```
#
# ``T^{\mathrm{out}}_{n\kappa}(r)``, which is related to the outer part of the orbital, is given in equation (283)
#
# ```math
# T^{\mathrm{out}}_{n\kappa}(r) = \beta_{n\kappa}
# + 2 \beta_{n\kappa} \frac{
#     a U(a+1, b+1, x_{\mathrm{out}}) + A \cdot (a+1) U(a+2, b+1, x_{\mathrm{out}})
# }{
#     U(a, b, x_{\mathrm{out}}) + A \cdot U(a + 1, b, x_{\mathrm{out}})
# }
# \\[1em]
# x_{\mathrm{out}} = 2 \beta_{n\kappa} r, \quad
# \beta_{n\kappa} = C^+ C^- / c, \quad
# C^{\pm} = \sqrt{c^2 \pm W_{n\kappa}}
# \\[1em]
# W_{n\kappa} = E + c^2, \quad
# v_{n\kappa} = \frac{1}{2} + \frac{Z W_{n\kappa}}{\beta_{n\kappa} c^2}, \quad
# A = \frac{Z}{\beta_{n\kappa}} + \kappa
# \\[1em]
# a = b/2 - v_{n\kappa}, \quad
# b = 2\gamma + 1, \quad
# \gamma = \sqrt{\kappa^2 - (Z/c)^2}
# ```
#
# We can get rid of the denominators, which may cause some unwanted singularities, and write down the following condition
#
# ```math
# \begin{aligned}
# f(E) & =
# \left(
#     \frac{\ell + 1 - \gamma}{R} + \beta_{n\kappa}
# \right) (U_{00} + A U_{10}) J_0 \\
# & - \alpha_{n\kappa} (U_{00} + A U_{10}) J_1 \\
# & + 2 \beta_{n\kappa} [a U_{11} + A (a + 1) U_{21}] ≡ 0
# \end{aligned}
# \\[1em]
# J_i = j_{\ell + 1}(x_{\mathrm{in}}), \quad
# U_{ij} = U(a + i, b + j, x_{\mathrm{out}})
# ```
#
# The special functions ``j_\ell(x)`` and ``U(a, b, x)`` are the spherical bessel and irregular confluent hypergeometric functions, respectively.
# They are evaluated using the corresponding functions from GSL ([`sf_bessel_jl`](https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=sf_bessel_jl#c.gsl_sf_bessel_jl), [`sf_hyperg_U`](https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=sf_hyperg_u#c.gsl_sf_hyperg_U)), via [GSL.jl](https://github.com/JuliaMath/GSL.jl).

## The following function can be used to visualize ``f``.
function plot_rel_fnc_f(; R, Z, κ, energies = nothing, plotfn=plot)
    es = if isnothing(energies)
        emin = hydrogenic_energy(DiracElectron, Z; n = 1, κ = κ)
        ## Excluding the first point, since at E = -Z/R, f is zero, which causes problems
        es = range(max(1.1emin, -Z/R), 0, length=502)[2:end-1]
    else
        energies
    end
    ys = (E -> AtomicMiscellany.f_andrae_dirac_topslice(Z=Z, R=R, κ=κ, E = E)).(es)
    any(iszero.(ys)) && @warn "Zeros in ys"
    plotfn(es, abs.(ys), yaxis=:log10, legend=false, c = (y -> y>0 ? 2 : 1).(ys))
    if isnothing(energies)
        absy_sorted = filter(!iszero, sort(abs.(ys)))
        ylims!(absy_sorted[1]/2, absy_sorted[490])
    end
    e0s = map(AtomicMiscellany.κ2ℓ(κ) .+ (1:6)) do n
        try
            hydrogenic_energy(DiracElectron, UniformShellNucleus(R), Z; n = n, κ = κ, verbose=false)
        catch e
            @warn "Unable to find root for n=$n" error=(e, catch_backtrace())
            NaN
        end
    end
    xlabel!("E"); ylabel!("log |f|")
    title!("Z=$(fmt(".2f", Z)), κ=$κ, R=$(fmt(".2e", R))")
    vline!(e0s, c=:gray)
    vline!([hydrogenic_energy(DiracElectron, Z; n = n, κ = κ) for n = 1:6], c=:black, ls=:dash)
    plot!(size=(800, 400))
end
nothing # hide

# The general behavior of ``f`` is very similar to the non-relativistic case.
# However, the expressions given in [`Andrae2000`](@ref) require us to evaluate ``\gamma = \sqrt{\kappa^2 - (Z/c)^2}``, so the expressions would become complex for ``Z > α^{-1}``. For that reason we restrict ourselves to ``Z \leq α^{-1}``.

plot(
    plot_rel_fnc_f(R = 1e-4, Z = 1, κ = -1),
    plot_rel_fnc_f(R = 5e-1, Z = 1, κ = -1),
    layout = (2, 1), size=(800, 600),
)

#-

plot_fnc_energies(DiracElectron, n = 1, κ = -1)

# As in the non-relativistic case, we run into issues with the lower energy zeros coming into the bracketing interval once the parameters become extreme enough (as we're limited in ``Z``, this primarily means ``R`` in this case).

plot_fnc_energies(DiracElectron, n = 2, κ = -1)
#-
plot_fnc_energies(DiracElectron, n = 3, κ = -1)
#-
plot_fnc_energies(DiracElectron, n = 2, κ = 1)

# We do run into trouble if we try to find the zeros for higher angular momenta:

plot_fnc_energies(DiracElectron, n = 2, κ = -2)
#-
plot_fnc_energies(DiracElectron, n = 3, κ = 2)

# This, turns out, is due to the limited resolution of 64-bit floats. We can zoom in on around the `1s` zero. We see that ``f(E)`` evaluates to the same value for many ``E`` values, implying a rounding problem. This, in turn, leads to the root-finder not being happy, since the bracketing interval does not encompass any positive ``f(E)`` values.

let R = 1e-3, Z = 1, κ = -2,
    e0 = hydrogenic_energy(DiracElectron, Z; n = 2, κ = κ),
    (emin, emax) = e0 .+ (-1, 1) .* 0.000_000_000_05
    plot_rel_fnc_f(R=R, Z=Z, κ=κ, energies = range(emin, emax, length=151), plotfn=scatter)
    xlims!(emin, emax)
end

# There are workarounds (modify the mathematical expressions for better numerical stability, use `BigFloat`s), but for now this is simply a limitation of the current implementation.
#
# Also, the non-relativistic version should also run into this issue at high enough ``\ell``.
