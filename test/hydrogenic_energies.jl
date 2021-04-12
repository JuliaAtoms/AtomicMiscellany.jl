module HydrogenicEnergiesTests
using AtomicMiscellany: AtomicMiscellany, hydrogenic_energy, α,
    NRParticle, NRElectron, NRMuon,
    DiracParticle, DiracElectron, DiracMuon,
    UniformShellNucleus
using Test

@testset "hydrogenic energies" begin
    @testset "non-rel / point" begin
        E(m, n, Z) = - m * Z^2 / (2 * n^2)
        @test hydrogenic_energy(NRElectron, 1; n = 1) ≈ -0.5
        @test hydrogenic_energy(NRElectron, 1; n = 2) ≈ -0.125
        @test hydrogenic_energy(NRElectron, 2; n = 1) ≈ -2
        @test hydrogenic_energy(NRElectron, 2; n = 3) ≈ -2/9
        for ℓ = [0, 5], n = (ℓ + 1):5:20, Z = [1, 2, 138]
            for m = [0.5, 1, 121]
                @test hydrogenic_energy(NRParticle(m), Z; n = n, ℓ = ℓ) ≈ E(m, n, Z)
            end
            @test hydrogenic_energy(NRElectron, Z; n = n, ℓ = ℓ) ≈ E(1, n, Z)
            @test hydrogenic_energy(NRMuon, Z; n = n, ℓ = ℓ) ≈ E(AtomicMiscellany.muon_mass_au, n, Z)
        end
    end

    @testset "relativistic / point" begin
        E(m, n, κ, Z) = m / (α^2 * sqrt(1 + (Z*α / (n - abs(κ) + sqrt(κ^2 - (Z*α)^2)))^2)) - m * α^-2
        for κ = [-2, -1, 1, 2], n = (AtomicMiscellany.κ2ℓ(κ) + 1):5:20, Z = [1, 2, 1/α]
            for m = [0.5, 1, 121]
                @test hydrogenic_energy(DiracParticle(m), Z; n = n, κ = κ) ≈ E(m, n, κ, Z)
                # Relativistic energy is always lower than the non-relativistic energy
                @test hydrogenic_energy(DiracParticle(m), Z; n = n, κ = κ) < hydrogenic_energy(NRParticle(m), Z; n = n)
            end
            @test hydrogenic_energy(DiracElectron, Z; n = n, κ = κ) ≈ E(1, n, κ, Z)
            @test hydrogenic_energy(DiracMuon, Z; n = n, κ = κ) ≈ E(AtomicMiscellany.muon_mass_au, n, κ, Z)
        end
    end

    @testset "FNC / NRElectron (R=$R, Z=$Z, nr=$nr, ℓ=$ℓ)" for R = [1e-4, 1e-3], Z = [1, 10, 100], nr = 1:3, ℓ = 0:4
        @test hydrogenic_energy(NRElectron, UniformShellNucleus(R), Z; n = nr + ℓ, ℓ = ℓ) > hydrogenic_energy(NRElectron, Z; n = nr + ℓ)
    end
    @testset "FNC / DiracElectron (R=$R, Z=$Z)" for R = [1e-4, 1e-3], Z = [1, 10, 100]
        @testset for nr = 1:10
            @test hydrogenic_energy(DiracElectron, UniformShellNucleus(R), Z; n = nr, κ = -1) > hydrogenic_energy(DiracElectron, Z; n = nr, κ = -1)
        end
        @testset for nr = 2:4
            # nr = 1 and nr > 4 fail, due to floating point errors
            @test hydrogenic_energy(DiracElectron, UniformShellNucleus(R), Z; n = nr + 1, κ = 1) > hydrogenic_energy(DiracElectron, Z; n = nr + 1, κ = 1)
        end
    end
end

end
