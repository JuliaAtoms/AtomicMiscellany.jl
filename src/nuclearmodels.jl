"""
    abstract type AbstractNuclearModel

Supertype for all (radially symmetric) nuclear models. The underlying charge distribution of
the nuclear model is assumed to be normalized to unity (i.e. corresponding to ``Z = 1``).

Each concrete nuclear model of type `T <: AbstractNuclearModel` should implement the
following method

* [`potential(::T, r::Real) -> Float64`](@ref potential):
  return the value of the potential at the radius `r`.

Additional, each nuclear model may also implement

* [`density(::T, r::Real) -> Float64`](@ref density):
  return the value of the (normalized) charge density at the radius `r`.

* [`rms(::T)`](@ref rms):
  return the root-mean-square radius ``\\sqrt{\\langle r^2 \\rangle}`` of the underlying
  charge distribution.

* [`from_rms(::Type{T}, rms)`](@ref from_rms):
  construct an instance of the nuclear model with the specified root-mean-square radius.

!!! note "Units"

    All quantities are assumed to be in atomic units.
"""
abstract type AbstractNuclearModel end
Base.Broadcast.broadcastable(x::AbstractNuclearModel) = Ref(x)

"""
    potential(::AbstractNuclearModel, r::Real) -> Float64

Return the value of the (normalized) charge density at the radius `r`.
"""
function potential end

"""
    density(::AbstractNuclearModel, r::Real) -> Float64

Return the value of the (normalized) charge density at the radius `r`.
"""
function density end

"""
    rms(::AbstractNuclearModel)

Return the root-mean-square radius ``\\sqrt{\\langle r^2 \\rangle}`` of the underlying
charge distribution.
"""
function rms end

"""
    from_rms(::Type{<:AbstractNuclearModel}, rms)

Construct an instance of the nuclear model with the specified root-mean-square radius.
"""
function from_rms end

"""
    struct PointNucleus <: AbstractNuclearModel

Represent a point nucleus with the potential

```math
V(r) = - Z / r,
\\quad
\\rho(r) = \\frac{Z}{4\\pi r^2} \\delta(r)
```

This nuclear model has no parameters. Also, this model does not define a method for
`density` due to the delta-function nature of the charge distribution.

# Constructors

* `PointNucleus()`: construct an instance of `PointNucleus`

# Examples

```jldoctest
julia> using AtomicMiscellany: PointNucleus, rms, potential

julia> rms(PointNucleus())
0.0

julia> potential(PointNucleus(), 2)
-0.5
```
"""
struct PointNucleus <: AbstractNuclearModel end
function potential(::PointNucleus, r::Real)
    r >= 0 || throw(DomainError(r, "r must be positive"))
    return -1/r
end
rms(::PointNucleus) = 0.0

"""
    struct UniformShellNucleus <: AbstractNuclearModel

Represents a nuclear model where all the charge is unformly distributed over an infinitely
thin shell at radius `R`. The nuclear potential and charge distribution are given by

```math
V(r) = \\begin{cases}
    -Z/R, & 0 \\leq r \\leq R \\\\
    -Z/r, & r > R
\\end{cases},
\\quad
\\rho(r) = \\frac{Z}{4\\pi r^2} \\delta(r - R)
```

Note that this model does not define a method for `density` due to the delta-function
nature of the charge distribution.

## Constructors

* `UniformShellNucleus(R::Real)`: construct a shell nucleus with radius `R`.

* `from_rms(UniformShellNucleus, rms)`: constructs a shell nucleus with the radius
  determined from the root-mean-square radius ``\\sqrt{\\langle r^2 \\rangle}`` of the
  charge distribution

  ```math
  R = \\sqrt{\\langle r^2 \\rangle}
  ```
"""
struct UniformShellNucleus <: AbstractNuclearModel
    R :: Float64
    function UniformShellNucleus(R::Real)
        R > 0 || throw(DomainError(R, "R must have a positive value"))
        new(R)
    end
end
from_rms(::Type{UniformShellNucleus}, rms) = UniformShellNucleus(rms)
function potential(m::UniformShellNucleus, r::Real)
    r >= 0 || throw(DomainError(r, "r must be positive"))
    r < m.R ? -1/m.R : -1/r
end
rms(m::UniformShellNucleus) = m.R

"""
    struct UniformSphericalNucleus <: AbstractNuclearModel

Represents a nuclear charge distribution where the charge is homogeneously distributed in
a sphere of radius ``R``. The potential and charge distributions are given by


```math
$(raw"""
V(r) = \begin{cases}
    -\frac{Z}{2R}\left[3 - \left(\frac{r}{R}\right)^2\right], & 0 \leq r \leq R \\
    -\frac{Z}{r}, & r > R
\end{cases},
\quad
\rho(r) = \begin{cases}
    \frac{3Z}{4\pi R^3}, & 0 \leq r \leq R \\
    0, & r > R
\end{cases}
""")
```

## Constructors

* `UniformSphericalNucleus(R::Real)`: construct a homogeneous spherical nucleus with radius
  `R`.

* `from_rms(UniformSphericalNucleus, rms)`: constructs a homogeneous spherical nucleus with
  the radius determined from the root-mean-square radius ``\\sqrt{\\langle r^2 \\rangle}``
  of the charge distribution

  ```math
  R = \\sqrt{\\frac{5}{3}} \\sqrt{\\langle r^2 \\rangle}
  ```

"""
struct UniformSphericalNucleus <: AbstractNuclearModel
    R :: Float64
    function UniformSphericalNucleus(R::Real)
        R > 0 || throw(DomainError(R, "R must have a positive value"))
        new(R)
    end
end
from_rms(::Type{UniformSphericalNucleus}, rms) = UniformSphericalNucleus(rms*sqrt(5/3))
function potential(m::UniformSphericalNucleus, r::Real)
    r >= 0 || throw(DomainError(r, "r must be positive"))
    R = m.R
    r < R ? ((r / R)^2 - 3) / (2 * R) : -1/r
end
function density(m::UniformSphericalNucleus, r::Real)
    r >= 0 || throw(DomainError(r, "r must be positive"))
    R = m.R
    r <= R ? 3 / (4π * R^3) : 0.0
end
rms(m::UniformSphericalNucleus) = m.R * sqrt(3/5)

"""
    struct GaussianNucleus <: AbstractNuclearModel

Represents a Gaussian-shaped nuclear charge distribution, size of which is determined by
the parameter ``R```. The potential and charge distributions are given by

```math
$(raw"""
V(r) = - \frac{Z}{r} \mathrm{erf}\left(\frac{r}{R}\right),
\quad
\rho(r) = \frac{Z}{\pi^{3/2} R^3}\exp\left[-\left(\frac{r}{R}\right)^2\right]
""")
```

## Constructors

* `GaussianNucleus(R::Real)`: constructs a Gaussian nuclear model of size `R`.

* `from_rms(UniformSphericalNucleus, rms)`: constructs a Gaussian nucleus with the ``R``
  parameter determined from the root-mean-square radius ``\\sqrt{\\langle r^2 \\rangle}``
  of the charge distribution

  ```math
  R = \\sqrt{\\frac{2}{3}} \\sqrt{\\langle r^2 \\rangle}
  ```

"""
struct GaussianNucleus <: AbstractNuclearModel
    R :: Float64

    function GaussianNucleus(R::Real)
        R > 0 || throw(DomainError(R, "R must have a positive value"))
        new(R)
    end
end
from_rms(::Type{GaussianNucleus}, rms) = GaussianNucleus(rms * sqrt(2/3))
function potential(m::GaussianNucleus, r::Real)
    r >= 0 || throw(DomainError(r, "r must be positive"))
    R = m.R
    # we need to special case r=0 due to the 1/r behaviour, even though the limit is fine.
    # We can get the correct value by considering the Taylor series of erf()
    # erf(x) = 2/√π (x - x^3 / 3 + ...)
    r > 0 ? - erf(r/R) / r : - 2/(R * sqrt(π))
end
function density(m::GaussianNucleus, r::Real)
    r >= 0 || throw(DomainError(r, "r must be positive"))
    R = m.R
    exp(-(r/R)^2) / (π^(3/2) * R^3)
end
rms(m::GaussianNucleus) = m.R * sqrt(3/2)

# """
#     struct FermiNucleus <: AbstractNuclearModel

# ...
# """
# struct FermiNucleus <: AbstractNuclearModel
#     c :: Float64
#     b :: Float64

#     function FermiNucleus() error("Not implemented.") end
# end
# from_rms(::Type{FermiNucleus}, rms) = PointNucleus() #FermiNucleus(0.0, 0.0)
# function potential(m::FermiNucleus, r::Real)
#     r >= 0 || throw(DomainError(r, "r must be positive"))
#     0.0
# end
# function density(m::FermiNucleus, r::Real)
#     r >= 0 || throw(DomainError(r, "r must be positive"))
#     0.0
# end
# rms(m::FermiNucleus) = m.R * sqrt(3/2)

"""
    struct Andrae2000

Type to dispatch on data from the following paper:

* Andrae, D. _"Finite Nuclear Charge Density Distributions in Electronic Structure
  Calculations for Atoms and Molecules."_ Physics Reports 336, no. 6 (October 2000):
  413–525. <https://doi.org/10.1016/S0370-1573(00)00007-7>.
"""
struct Andrae2000 end
