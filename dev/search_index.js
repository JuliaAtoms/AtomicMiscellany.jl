var documenterSearchIndex = {"docs":
[{"location":"implementation/nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"CurrentModule = AtomicMiscellany","category":"page"},{"location":"implementation/nuclearmodels/#Nuclear-models","page":"Nuclear models","title":"Nuclear models","text":"","category":"section"},{"location":"implementation/nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"Notes on the internal implementation details of the nuclear models. Anything discussed here is not considered to be part of the public API.","category":"page"},{"location":"implementation/nuclearmodels/#Density-vs.-potential","page":"Nuclear models","title":"Density vs. potential","text":"","category":"section"},{"location":"implementation/nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"The nuclear charge density and the corresponding potential are related to each other via the Poisson equation. In radial coordinates, for a spherically symmetric density, the potential at r, as per equation (22) in Andrae2000, is given by the following integral","category":"page"},{"location":"implementation/nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"V(r) = - frac4pirleft\n    int_0^r s^2 rho(s) ds\n    + r int_r^infty s rho(s) ds\nright","category":"page"},{"location":"implementation/nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"In the tests, we use this to numerically verify that the implementations for density and potential actually match, as they are generally implemented independently based on analytical expressions.","category":"page"},{"location":"implementation/nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"The numerical integration can be done using the following routine","category":"page"},{"location":"implementation/nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"using QuadGK\nfunction potential_quadk(ρ, r; rmax = Inf, rscale = 1.0)\n    # Based on equation (22) in Andrae. The potential is split into two separate integrals\n    # on the intervals [0, r] and [r, ∞).\n    #\n    # However, often the density can have its non-trivial parts at very small r values,\n    # which means that the adaptive quadgk algorithm might not realize that it even has\n    # non-zero parts. To work around that, we do 2 things:\n    #\n    # 1. We scale the r-axis with rscale (so r = rscale*t), and actually integrate in t.\n    #    In practice, we'll use the RMS value for the scaling, since we'd expect that to\n    #    correspond to the scale at which the density has a non-trivial shape.\n    # 2. If r is much greater than rscale, we split the [0, r] integration into two\n    #    sub-intervals, hoping that we capture all the non-trivial part in the first\n    #    interval.\n    #\n    # The change of variable r = rscale*t leads to\n    #\n    #   ∫ f(r) dr = rscale ∫ f(rscale*t) dt\n    #\n    # So we have to also multiply the resulting integral with rscale.\n    _I1(t) = rscale^3 * t^2 * ρ(rscale*t)\n    I1, δ1 = min(r, rmax)/rscale > 10.0 ? quadgk(_I1, 0, 10.0, min(r, rmax)/rscale) :\n        quadgk(_I1, 0, min(r, rmax)/rscale)\n    I2, δ2 = (r < rmax) ? quadgk(t -> rscale^2 * t * ρ(rscale*t), r/rscale, rmax/rscale) : (0.0, 0.0)\n    # add up the subintegral based on equation (22):\n    -4π*(I1/r + I2)\nend\nnothing # hide","category":"page"},{"location":"implementation/nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"As an example, let's plot a few arbitrary cases","category":"page"},{"location":"implementation/nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"using AtomicMiscellany.NuclearModels\nusing Plots: plot, plot!\nrs = range(0, 1, length=101)\nplot(legend=:bottomright)\nlet nm = UniformSphericalNucleus(0.6)\n    plot!(\n        rs, potential_quadk.(r -> density(nm, r), rs, rmax = 0.6, rscale = rms(nm)),\n        label = \"Uniform sphere\", c = 1, lw = 3,\n    )\n    plot!(rs, (r -> potential(nm, r)).(rs), c=:black, ls=:dash, label=nothing)\nend\nlet nm = GaussianNucleus(0.25)\n    plot!(\n        rs, potential_quadk.(r -> density(nm, r), rs, rscale = rms(nm)),\n        label = \"Gaussian\", c = 2, lw = 3,\n    )\n    plot!(rs, (r -> potential(nm, r)).(rs), c=:black, ls=:dash, label=\"Analytic\")\nend","category":"page"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"CurrentModule = AtomicMiscellany","category":"page"},{"location":"nuclearmodels/#Nuclear-models","page":"Nuclear models","title":"Nuclear models","text":"","category":"section"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"Each nuclear model is implemented as a subtype of AbstractNuclearModel. An instance of each nuclear model is assumed to be normalized to Z = 1.","category":"page"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"AbstractNuclearModel","category":"page"},{"location":"nuclearmodels/#AtomicMiscellany.AbstractNuclearModel","page":"Nuclear models","title":"AtomicMiscellany.AbstractNuclearModel","text":"abstract type AbstractNuclearModel\n\nSupertype for all (radially symmetric) nuclear models. The underlying charge distribution of the nuclear model is assumed to be normalized to unity (i.e. corresponding to Z = 1).\n\nEach concrete nuclear model of type T <: AbstractNuclearModel should implement the following method\n\npotential(::T, r::Real) -> Float64: return the value of the potential at the radius r.\n\nAdditional, each nuclear model may also implement\n\ndensity(::T, r::Real) -> Float64: return the value of the (normalized) charge density at the radius r.\nrms(::T): return the root-mean-square radius sqrtlangle r^2 rangle of the underlying charge distribution.\nfrom_rms(::Type{T}, rms): construct an instance of the nuclear model with the specified root-mean-square radius.\n\nnote: Units\nAll quantities are assumed to be in atomic units.\n\n\n\n\n\n","category":"type"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"The package implements various nuclear models, illustrated by the following plot. The root-mean-square radius of the charge density is the same for each mode to keep the potentials comparable.","category":"page"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"using AtomicMiscellany.NuclearModels, Plots\nrs = 10 .^ range(-6, -2, length=501)\nplot(\n    title = \"Nuclear potentials with RMS = 10⁻⁴ a.u.\",\n    xlabel = \"r (a.u.)\", ylabel = \"Nuclear potential V\",\n    xaxis=:log10, ylims = (-2.1e4, 0.1e4), legend=:bottomright,\n)\nV(M) = r -> potential(M, r)\nplot!(rs, V(PointNucleus()).(rs), label = \"PointNucleus\")\nplot!(rs, V(from_rms(UniformShellNucleus, 1e-4)).(rs), label = \"UniformShellNucleus\")\nplot!(rs, V(from_rms(UniformSphericalNucleus, 1e-4)).(rs), label = \"UniformSphericalNucleus\")\nplot!(rs, V(from_rms(GaussianNucleus, 1e-4)).(rs), label = \"GaussianNucleus\")\n#plot!(rs, V(from_rms(FermiNucleus, 1e-4)).(rs), label = \"FermiNucleus\") # hide\nplot!() # hide","category":"page"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"PointNucleus\nUniformShellNucleus\nUniformSphericalNucleus\nGaussianNucleus","category":"page"},{"location":"nuclearmodels/#AtomicMiscellany.PointNucleus","page":"Nuclear models","title":"AtomicMiscellany.PointNucleus","text":"struct PointNucleus <: AbstractNuclearModel\n\nRepresent a point nucleus with the potential\n\nV(r) = - Z  r\nquad\nrho(r) = fracZ4pi r^2 delta(r)\n\nThis nuclear model has no parameters. Also, this model does not define a method for density due to the delta-function nature of the charge distribution.\n\nConstructors\n\nPointNucleus(): construct an instance of PointNucleus\n\nExamples\n\njulia> using AtomicMiscellany: PointNucleus, rms, potential\n\njulia> rms(PointNucleus())\n0.0\n\njulia> potential(PointNucleus(), 2)\n-0.5\n\n\n\n\n\n","category":"type"},{"location":"nuclearmodels/#AtomicMiscellany.UniformShellNucleus","page":"Nuclear models","title":"AtomicMiscellany.UniformShellNucleus","text":"struct UniformShellNucleus <: AbstractNuclearModel\n\nRepresents a nuclear model where all the charge is unformly distributed over an infinitely thin shell at radius R. The nuclear potential and charge distribution are given by\n\nV(r) = begincases\n    -ZR  0 leq r leq R \n    -Zr  r  R\nendcases\nquad\nrho(r) = fracZ4pi r^2 delta(r - R)\n\nNote that this model does not define a method for density due to the delta-function nature of the charge distribution.\n\nConstructors\n\nUniformShellNucleus(R::Real): construct a shell nucleus with radius R.\nfrom_rms(UniformShellNucleus, rms): constructs a shell nucleus with the radius determined from the root-mean-square radius sqrtlangle r^2 rangle of the charge distribution\nR = sqrtlangle r^2 rangle\n\n\n\n\n\n","category":"type"},{"location":"nuclearmodels/#AtomicMiscellany.UniformSphericalNucleus","page":"Nuclear models","title":"AtomicMiscellany.UniformSphericalNucleus","text":"struct UniformSphericalNucleus <: AbstractNuclearModel\n\nRepresents a nuclear charge distribution where the charge is homogeneously distributed in a sphere of radius R. The potential and charge distributions are given by\n\nV(r) = begincases\n    -fracZ2Rleft3 - left(fracrRright)^2right  0 leq r leq R \n    -fracZr  r  R\nendcases\nquad\nrho(r) = begincases\n    frac3Z4pi R^3  0 leq r leq R \n    0  r  R\nendcases\n\n\nConstructors\n\nUniformSphericalNucleus(R::Real): construct a homogeneous spherical nucleus with radius R.\nfrom_rms(UniformSphericalNucleus, rms): constructs a homogeneous spherical nucleus with the radius determined from the root-mean-square radius sqrtlangle r^2 rangle of the charge distribution\nR = sqrtfrac53 sqrtlangle r^2 rangle\n\n\n\n\n\n","category":"type"},{"location":"nuclearmodels/#AtomicMiscellany.GaussianNucleus","page":"Nuclear models","title":"AtomicMiscellany.GaussianNucleus","text":"struct GaussianNucleus <: AbstractNuclearModel\n\nRepresents a Gaussian-shaped nuclear charge distribution, size of which is determined by the parameter R`. The potential and charge distributions are given by\n\nV(r) = - fracZr mathrmerfleft(fracrRright)\nquad\nrho(r) = fracZpi^32 R^3expleft-left(fracrRright)^2right\n\n\nConstructors\n\nGaussianNucleus(R::Real): constructs a Gaussian nuclear model of size R.\nfrom_rms(UniformSphericalNucleus, rms): constructs a Gaussian nucleus with the R parameter determined from the root-mean-square radius sqrtlangle r^2 rangle of the charge distribution\nR = sqrtfrac23 sqrtlangle r^2 rangle\n\n\n\n\n\n","category":"type"},{"location":"nuclearmodels/#Root-mean-square-radius","page":"Nuclear models","title":"Root-mean-square radius","text":"","category":"section"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"The following methods are for determining the RMS radius of nuclei, loosely based on experimental data.","category":"page"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"rms(::Type{JohnsonSoff1985}, ::Real)","category":"page"},{"location":"nuclearmodels/#AtomicMiscellany.rms-Tuple{Type{AtomicMiscellany.JohnsonSoff1985}, Real}","page":"Nuclear models","title":"AtomicMiscellany.rms","text":"rms(JohnsonSoff1985, A::Real) -> Float64\n\nReturns the root-mean-square radius (in atomic units) of a nucleus based on the fit from the JohnsonSoff1985 paper\n\nrmrms(A) = (0836 A^13 + 0570)rmfm\n\nwhere A is the atomic mass number of the isotope.\n\njulia> using AtomicMiscellany: rms, JohnsonSoff1985\n\njulia> rms(JohnsonSoff1985, 238)\n0.00010867476884785438\n\nnote: Accuracy and validity\nThe paper states that the fit is for A  9, and it appears that the largest A value in the data used for the fit was about 250. The implementation, however, will work with any positive A value.The accuracy is stated to be 005rmfm i.e. 95 times 10^-7rmau.\n\n\n\n\n\n","category":"method"},{"location":"nuclearmodels/#Utilities","page":"Nuclear models","title":"Utilities","text":"","category":"section"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"potential\ndensity\nrms\nfrom_rms\nAtomicMiscellany.NuclearModels","category":"page"},{"location":"nuclearmodels/#AtomicMiscellany.potential","page":"Nuclear models","title":"AtomicMiscellany.potential","text":"potential(::AbstractNuclearModel, r::Real) -> Float64\n\nReturn the value of the (normalized) charge density at the radius r.\n\n\n\n\n\n","category":"function"},{"location":"nuclearmodels/#AtomicMiscellany.density","page":"Nuclear models","title":"AtomicMiscellany.density","text":"density(::AbstractNuclearModel, r::Real) -> Float64\n\nReturn the value of the (normalized) charge density at the radius r.\n\n\n\n\n\n","category":"function"},{"location":"nuclearmodels/#AtomicMiscellany.rms","page":"Nuclear models","title":"AtomicMiscellany.rms","text":"rms(::AbstractNuclearModel)\n\nReturn the root-mean-square radius sqrtlangle r^2 rangle of the underlying charge distribution.\n\n\n\n\n\n","category":"function"},{"location":"nuclearmodels/#AtomicMiscellany.from_rms","page":"Nuclear models","title":"AtomicMiscellany.from_rms","text":"from_rms(::Type{<:AbstractNuclearModel}, rms)\n\nConstruct an instance of the nuclear model with the specified root-mean-square radius.\n\n\n\n\n\n","category":"function"},{"location":"nuclearmodels/#AtomicMiscellany.NuclearModels","page":"Nuclear models","title":"AtomicMiscellany.NuclearModels","text":"Can be used to bring various bindings related to nuclear models into the namespace.\n\njulia> using AtomicMiscellany.NuclearModels\n\njulia> parentmodule(rms), parentmodule(AbstractNuclearModel)\n(AtomicMiscellany, AtomicMiscellany)\n\n\n\n\n\n","category":"module"},{"location":"nuclearmodels/#References","page":"Nuclear models","title":"References","text":"","category":"section"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"Andrae2000\nJohnsonSoff1985","category":"page"},{"location":"nuclearmodels/#AtomicMiscellany.Andrae2000","page":"Nuclear models","title":"AtomicMiscellany.Andrae2000","text":"struct Andrae2000\n\nType to dispatch on data from the following paper:\n\nAndrae, D. \"Finite Nuclear Charge Density Distributions in Electronic Structure Calculations for Atoms and Molecules.\" Physics Reports 336, no. 6 (October 2000): 413–525. https://doi.org/10.1016/S0370-1573(00)00007-7.\n\n\n\n\n\n","category":"type"},{"location":"nuclearmodels/#AtomicMiscellany.JohnsonSoff1985","page":"Nuclear models","title":"AtomicMiscellany.JohnsonSoff1985","text":"struct JohnsonSoff1985\n\nType to dispatch on data from the following paper:\n\nW.R. Johnson, Gerhard Soff, \"The lamb shift in hydrogen-like atoms, 1 ⩽ Z ⩽ 110\", Atomic Data and Nuclear Data Tables, Volume 33, Issue 3, 1985, Pages 405-446, https://doi.org/10.1016/0092-640X(85)90010-5\n\nSee also: rms(::Type{JohnsonSoff1985}, ::Real).\n\n\n\n\n\n","category":"type"},{"location":"nuclearmodels/#Index","page":"Nuclear models","title":"Index","text":"","category":"section"},{"location":"nuclearmodels/","page":"Nuclear models","title":"Nuclear models","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = AtomicMiscellany","category":"page"},{"location":"#AtomicMiscellany","page":"Home","title":"AtomicMiscellany","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"AtomicMiscellany","category":"page"},{"location":"#AtomicMiscellany.AtomicMiscellany","page":"Home","title":"AtomicMiscellany.AtomicMiscellany","text":"A gallimaufry of disparate utilities useful in atomic physics that are too small and/or simple to deserve a package of their own.\n\nBy default, AtomicMiscellany does not export anything (i.e. doing using AtomicMiscellany does not introduce any new names into the namespace). This is because the package contains a lot of unrelated bindings and would is most cases pollute the namespace. Instead, you can either being the necessary binding into the namespace explicitly or fully qualify them at the use-site with the package name.\n\nThe package also provides submodules to bring related subsets bindings into namespace:\n\nAtomicMiscellany.NuclearModels\n\n\n\n\n\n","category":"module"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}