"""
    abstract type AbstractParticle

Subtypes of `AbstractParticle` represent various particles and/or Hamiltonians, such as
particles described the Dirac equation or the non-relativistic Schrödinger equation.

The subtypes are used to determine which physical particle species is meant when dispatching
on generic functions calculating some properties of or related to the particle.
"""
abstract type AbstractParticle end
Base.Broadcast.broadcastable(x::AbstractParticle) = Ref(x)

"""
    struct DiracParticle <: AbstractParticle

Represents a relativistic particle, the free-particle dynamics of which are described by the
Dirac equation.

```math
(\\gamma^\\mu \\partial_\\mu - m) \\Psi = 0
```

Note that we set the zero of the energy to ``+mc^2``.
"""
struct DiracParticle <: AbstractParticle
    mass :: Float64
    function DiracParticle(mass)
        mass > 0 || throw(DomainError(mass, "mass must be positive"))
        new(mass)
    end
end
"""
    const DiracElectron :: DiracParticle

A relativistic Dirac fermion ([`DiracParticle`](@ref)) with the mass of an electron.
"""
const DiracElectron = DiracParticle(1.0)
"""
    const DiracMuon :: DiracParticle

A relativistic Dirac fermion ([`DiracParticle`](@ref)) with the mass of a muon.
"""
const DiracMuon = DiracParticle(muon_mass_au)

"""
    struct NRParticle <: AbstractParticle

Represents a non-relativistic particle, described by the non-relativistic Schrödinger
equation.

```math
i \\frac{\\partial\\Psi}{\\partial t} = - \\frac{\\nabla^2}{2 m^2} \\Psi
```
"""
struct NRParticle <: AbstractParticle
    mass :: Float64
    function NRParticle(mass)
        mass > 0 || throw(DomainError(mass, "mass must be positive"))
        new(mass)
    end
end
"""
    const NRElectron

A non-relativistic particle ([`NRParticle`](@ref)) with the mass of an electron.
"""
const NRElectron = NRParticle(1.0)
"""
    const NRMuon

A non-relativistic particle ([`NRParticle`](@ref)) with the mass of a muon.
"""
const NRMuon = NRParticle(muon_mass_au)

"""
    hydrogenic_energy(
        dynamics::AbstractParticle[, nucleus::AbstractNuclearModel], Z;
        qnumbers...
    ) -> Float64

Returns the energy (in atomic units) of a hydrogenic bound state for given dynamics
(particle), central potential generated by a nuclear charge distribution, and quantum numbers.

Each nuclear model is, by the definition of [`AbstractNuclearModel`](@ref), normalized to unity.
This means that the nuclear charge ``Z`` needs to be passed as a separate argument.
Also, the second argument (`nucleus`) can be omitted and it defaults to [`PointNucleus`](@ref).

## Quantum numbers

The keyword arguments (`qnumbers...`) can be used to specify all the necessarily quantum
numbers to uniquely identify the state. In some cases, some quantum numbers may be omitted.

In the non-relativistic case we need to specify `n` and `ℓ`, whereas in the relativistic case
it is `n` and `κ`.
"""
function hydrogenic_energy end
# We let the second argument default to `PointNucleus`:
hydrogenic_energy(p::AbstractParticle, Z::Number; kwargs...) = hydrogenic_energy(p, PointNucleus(), Z; kwargs...)

function hydrogenic_energy(p::NRParticle, ::PointNucleus, Z::Real; n::Integer, ℓ = 0, m = 0, ms = 0)
    n >= 1 || throw(DomainError(n, "n must be a positive integer"))
    - p.mass * Z^2 / (2 * n^2)
end

function hydrogenic_energy(p::DiracParticle, ::PointNucleus, Z::Number; n::Integer, κ::Integer, mj = 0)
    n >= 1 || throw(DomainError(n, "n must be a positive integer"))
    if Z isa Real
        Z < 0 && throw(DomainError(Z, "0 ≤ Z ≤ α⁻¹ for real Z"))
        Z > 1/α && throw(DomainError(Z, "0 ≤ Z ≤ α⁻¹ for real Z, pass as complex(Z) to evaluate at Z > α⁻¹"))
    end
    # Note: using this expression is not as numerically precise as it could probably be, in
    # terms of floating point precision. It leads to an error of about 1e-12.
    (1 / sqrt(1 + (Z*α / (n - abs(κ) + sqrt(κ^2 - (Z*α)^2)))^2) - 1) * p.mass * α^-2
end

function hydrogenic_energy(p::NRParticle, nm::UniformShellNucleus, Z::Real; n::Integer, ℓ::Integer, m = 0, ms = 0, verbose=false)
    n >= 1 || throw(DomainError(n, "n must be a positive integer"))
    f = E -> f_andrae_nonrel_topslice(; Z=Z, R=nm.R, r=nm.R, ℓ=ℓ, E = E)
    # Pick the initial energy range to be between the point energies for n and n+1:
    emin, emax = hydrogenic_energy(p, Z; n = n), hydrogenic_energy(p, Z; n = n + 1)
    # In the expressions, we need to evaluate sqrt(2*(E+Z/R)), so E ≥ -Z/R for f to even
    # make sense. So we must constrain the bounding box when Z becomes large enough.
    emin = max(emin, -Z/nm.R)
    find_zero(f, (emin, emax), verbose = verbose)
end

# See docs/src/implementation/hydrogenic.jl for documentation.
function f_andrae_nonrel_topslice(; Z, R, ℓ, r, E)
    # Inner variables
    α = sqrt(2 * (E + Z/R))
    xin = α*r
    J1 = sf_bessel_jl(ℓ + 1, xin)
    J2 = sf_bessel_jl(ℓ, xin)
    # Outer variables
    βnℓ = sqrt(-2*E)
    νnℓ = Z/βnℓ
    b = 2*(ℓ+1)
    a = b/2 - νnℓ
    xout = 2*βnℓ*r
    U1 = try
        sf_hyperg_U(a+1, b+1, xout)
    catch e
        @error "sf_hyperg_U(a+1, b+1, xout)" a b xout r
        rethrow(e)
    end
    U2 = try
        sf_hyperg_U(a, b, xout)
    catch e
        @error "sf_hyperg_U(a, b, xout)" a b xout r
        rethrow(e)
    end
    # Multiply both sides and subtract
    βnℓ*J2*U2 + 2*βnℓ*a*J2*U1 - α*J1*U2
end

function hydrogenic_energy(p::DiracParticle, nm::UniformShellNucleus, Z::Real; n::Integer, κ::Integer, mj = 0, verbose=false)
    n >= 1 || throw(DomainError(n, "n must be a positive integer"))
    # We'll do root finding in a bracketing box, so we first need to decide the bracketing
    # interval. As we know that the energies always slightly shift up, we choose the
    # interval to be between the n and n+1 point energies. In most cases, this will
    # guarantee that we will have exactly one zero within the interval. Where this will fail
    # is when the FNC shift is so large that energy is greater than the next PNC energy.
    emin = hydrogenic_energy(p, Z; κ = κ, n = n)
    emax = hydrogenic_energy(p, Z; κ = κ, n = n + 1)
    # We shift the bracketing interval a bit because of the numerical errors in
    # hydrogenic_energy. With 64-bit floats, it's accuracy is only about 1e-12, so we are
    # conservative here and shift everything by 1e-10. This is still many orders of magnitude
    # smaller than the energy difference between any two consective n values (for low enough
    # n values), so this is pretty safe thing to do.
    #
    # However, the result is that we are guaranteed to have a sane bracketing interval,
    # whereas without this shift we sometimes place the start or end of the interval _after_
    # the corresponding zero, as opposed to having it before it.
    emin -= 1e-10
    emax -= 1e-10
    # In the expressions, we need to evaluate α_nκ, which require E ≥ -Z/R for f to not be
    # complex. So we constrain the bounding box when Z becomes large enough.
    emin = max(emin, -Z/nm.R)
    f = E -> f_andrae_dirac_topslice(; Z = Z, R = nm.R, κ = κ, E = E)
    find_zero(f, (emin, emax), verbose = verbose)
end

# See docs/src/implementation/hydrogenic.jl for documentation.
function f_andrae_dirac_topslice(; Z, R, κ, E)
    # Note: E = 0 is defined to be mc², so the energy will always be negative for bound
    # states
    c = α^-1 # since we are in atomic units
    ℓ = κ2ℓ(κ)
    # Following are from equation (282)
    #
    # Wnκ = E + c^2
    # C⁺, C⁻ = sqrt(c^2 + Wnκ), sqrt(c^2 - Wnκ)
    # βnκ = C⁺ * C⁻ / c
    #
    # However, as the expressions in the paper are subject to numerical roundoff errors that
    # limit the resolution of E to about 1e-12 (for 64-bit floats), we slightly rearrange
    # the equations for βnκ and also νnκ below.
    Ec⁻² = E/c^2
    βnκ = sqrt(-E * (2 + Ec⁻²))
    x_out = 2*βnκ*R # x at r = R for outer part
    A = Z/βnκ + κ
    γ = sqrt(κ^2 - (Z/c)^2)
    b = 2γ +  1
    # We also modify νnκ
    #
    # νnκ = 0.5 + Z*Wnκ/(βnκ * c^2)
    #
    # to reduce roundoff errors.
    νnκ = 0.5 + (Z/βnκ) * (Ec⁻² + 1)
    a = b/2 - νnκ
    # From 283
    U00 = sf_hyperg_U(a, b, x_out)
    U10 = sf_hyperg_U(a + 1, b, x_out)
    U11 = sf_hyperg_U(a + 1, b + 1, x_out)
    U21 = sf_hyperg_U(a + 2, b + 1, x_out)
    Ud = U00 + A*U10 # denomimator in (283)
    # From 273 & 274
    αnκ = sqrt(2*(E + Z/R) + ((E + Z/R)/c)^2)
    x_in = αnκ*R # x at r = R for inner part
    J0 = sf_bessel_jl(ℓ, x_in)
    J1 = sf_bessel_jl(ℓ + 1, x_in)
    # Combining all of them now:
    ((ℓ + 1 - γ)/R + βnκ)*Ud*J0 - αnκ*J1*Ud + 2*βnκ*(a*U11 + A*(a + 1)*U21)*J0
end

# From AtomicLevels.jl (MIT license)
function κ2ℓ(κ::Integer)
    κ == zero(κ) && throw(ArgumentError("κ can not be zero"))
    (κ < 0) ? -(κ + 1) : κ
end
