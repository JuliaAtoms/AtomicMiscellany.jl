"""
    struct JohnsonSoff1985

Type to dispatch on data from the following paper:

* W.R. Johnson, Gerhard Soff, _"The lamb shift in hydrogen-like atoms, 1 ⩽ Z ⩽ 110"_, Atomic
  Data and Nuclear Data Tables, Volume 33, Issue 3, 1985, Pages 405-446, <https://doi.org/10.1016/0092-640X(85)90010-5>

See also: [`rms(::Type{JohnsonSoff1985}, ::Real)`](@ref).
"""
struct JohnsonSoff1985 end

"""
    rms(JohnsonSoff1985, A::Real) -> Float64

Returns the root-mean-square radius (in atomic units) of a nucleus based on the fit from the
[`JohnsonSoff1985`](@ref) paper

```math
\\rm{rms}(A) = (0.836 A^{1/3} + 0.570)~\\rm{fm}
```

where `A` is the atomic mass number of the isotope.

```jldoctest
julia> using AtomicMiscellany: rms, JohnsonSoff1985

julia> rms(JohnsonSoff1985, 238)
0.00010867476884785438
```

!!! note "Accuracy and validity"

    The paper states that the fit is for ``A > 9``, and it appears that the largest ``A``
    value in the data used for the fit was about ``250``. The implementation, however, will
    work with any positive ``A`` value.

    The accuracy is stated to be ``0.05~\\rm{fm}`` i.e. ``9.5 \\times 10^{-7}~\\rm{a.u.}``.
"""
function rms(::Type{JohnsonSoff1985}, A::Real)
    A >= 0 || throw(DomainError(A, "Atomic mass number A must be strictly positive"))
    rms_fm = 0.836 * A^(1/3) + 0.570 # equation (20) from the paper
    # convert to atomic units
    rms_fm * (1e-15 / bohr_in_m)
end

# According to CODATA 2018: 5.29177210903(80)×10−11
const bohr_in_m = 5.29177210903e-11
