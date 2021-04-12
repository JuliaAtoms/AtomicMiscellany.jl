# Extra constants, expanding PhysicalConstants
@constant(
    MuonMass, m_μ, "Muon mass",
    1.883_531_627e-28,
    BigFloat(1_883_531_627)/BigFloat(10_000_000_000_000_000_000_000_000_000_000_000_000),
    kg,
    4.2e-36,
    BigFloat(42)/BigFloat(10_000_000_000_000_000_000_000_000_000_000_000_000),
    "CODATA 2018"
)

const muon_mass_au = AtomicMiscellany.MuonMass / CODATA2018.ElectronMass

const α = real(CODATA2018.FineStructureConstant)
