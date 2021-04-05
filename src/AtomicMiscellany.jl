"""
A gallimaufry of disparate utilities useful in atomic physics that are too small and/or
simple to deserve a package of their own.

By default, [`AtomicMiscellany`](@ref) does not export anything (i.e. doing
`using AtomicMiscellany` does not introduce any new names into the namespace). This is
because the package contains a lot of unrelated bindings and would is most cases pollute the
namespace. Instead, you can either being the necessary binding into the namespace explicitly
or fully qualify them at the use-site with the package name.

The package also provides submodules to bring related subsets bindings into namespace:

* [`AtomicMiscellany.NuclearModels`](@ref)
"""
module AtomicMiscellany
using SpecialFunctions: erf

include("nuclearmodels.jl")
include("johnsonsoff1985.jl")

"""
Can be used to bring various bindings related to nuclear models into the namespace.

```jldoctest
julia> using AtomicMiscellany.NuclearModels

julia> parentmodule(rms), parentmodule(AbstractNuclearModel)
(AtomicMiscellany, AtomicMiscellany)
```
"""
module NuclearModels
    import AtomicMiscellany: potential, density, rms, from_rms, AbstractNuclearModel,
        PointNucleus, UniformShellNucleus, UniformSphericalNucleus, GaussianNucleus,
        JohnsonSoff1985
    export potential, density, rms, from_rms, AbstractNuclearModel,
        PointNucleus, UniformShellNucleus, UniformSphericalNucleus, GaussianNucleus,
        JohnsonSoff1985
end

end
