module NuclearModelTests
using AtomicMiscellany.NuclearModels
using Test
using QuadGK

maxnonzero(::AbstractNuclearModel) = Inf
maxnonzero(m::UniformSphericalNucleus) = m.R

@testset "Nuclear models" begin
    # Test that any AbstractNuclearModel broadcasts like a scalar
    @test length(potential.(PointNucleus(), [1,2,3])) == 3
    @test length(density.(GaussianNucleus(1), [1,2,3,4])) == 4

    # Check that if we construct the nuclear models using the from_rms method, that the
    # resulting objects have the correct RMS value.
    function check_rms_constructor(T::Type{<:AbstractNuclearModel})
        @testset "check_rms_constructor: $T" begin
            for _rms in [3.1e-6, 1e-4, 0.4, 1.0, 10.0]
                m = from_rms(T, _rms)
                @test rms(m) ≈ _rms
            end
        end
    end

    @testset "from_rms" begin
        @test_throws MethodError from_rms(PointNucleus, 0.0)
        check_rms_constructor(UniformShellNucleus)
        check_rms_constructor(UniformSphericalNucleus)
        check_rms_constructor(GaussianNucleus)
        #check_rms_constructor(FermiNucleus)
    end

    # Use direct integration to calculate the RMS values and compare against the implemented
    # equations.
    function compare_integrals(m::AbstractNuclearModel)
        rmax, rscale = maxnonzero(m), rms(m)
        @testset "compare_integrals: $m" begin
            # We need to scale and split the integrals: see the comments for potential_quadk()
            # below.
            ρ₀, δ = let f(t) = 4π * rscale^3 * t^2 * density(m, rscale*t)
                rmax/rscale < 10.0 ? quadgk(f, 0, rmax/rscale) : quadgk(f, 0, 10.0, rmax/rscale)
            end
            @test ρ₀ ≈ 1.0

            r², δ = let f(t) = 4π * rscale^5 * t^4 * density(m, rscale*t)
                rmax/rscale < 10.0 ? quadgk(f, 0, rmax/rscale) : quadgk(f, 0, 10.0, rmax/rscale)
            end
            @debug "$r² ± $δ" sqrt(r²) rms(m) rms(m)^2
            @test sqrt(r²) ≈ rms(m)
        end
    end

    @testset "r^k integrals" begin
        for R in [1, 1e-4, 5.1e-6, 35, π, 5π^4]
            compare_integrals(UniformSphericalNucleus(R))
            compare_integrals(GaussianNucleus(R))
            #compare_integrals(FermiNucleus(R))
        end
    end

    # Test to make sure that the density and the potential correspond to each other, for the
    # nuclear models that we have the density implemented. We'll use QuadGK to numerically
    # integrate the potential values from the density directly and then check against the
    # implementation of the potential.
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

    test_potential(m::AbstractNuclearModel) = @testset "test_potential: $m" begin
        for r0 in 10 .^ range(-6, 2, length=21)
            V = potential_quadk(r -> density(m, r), r0, rmax = maxnonzero(m), rscale = rms(m))
            @test V ≈ potential(m, r0)
        end
    end

    @testset "V(r) from ρ(r)" begin
        for R in [1, 1e-4, 5.1e-6, 35, π, 5π^4]
            test_potential(UniformSphericalNucleus(R))
            test_potential(GaussianNucleus(R))
            #test_potential(FermiNucleus(R))
        end
    end
end

end
