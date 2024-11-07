using Test
using Unitful
using JLD2
using Statistics
using CanonicalConfigurations

# Upload test data that looks good
# check velocities in classical come from boltzmann
# check amplitudes approach eachother for high T

test_freqs_sq, test_phi, test_dynmat = 
    load(joinpath(@__DIR__, "SW_1300K_TestData.jld2"), "freqs_sq", "phi", "dynmat")
test_kB = uconvert(u"eV / K", Unitful.k)
test_hbar = uconvert(u"eV*s", Unitful.ħ)
test_freq_units = sqrt(u"eV / Å^2 / u")

test_freqs = Float32.(sqrt.(test_freqs_sq)) * test_freq_units
test_temp = 1300u"K"

N_atoms_test = 216
test_masses = Float32.(28.085u"u" .* ones(N_atoms_test))
extended_test_masses = collect(Iterators.flatten(Iterators.repeated(el, 3) for el in test_masses))

function KE(v, m)
    return 0.5 * sum(v.^2 .* m)
end

@testset "ClassicalKineticEnergy" begin
    # Test that average kinetic energy of each DoF is 0.5*kB*T
    n_configs = 250_000
    settings = ClassicalConfigSettings(test_kB, n_configs, test_temp)

    velos = canonical_velocities(settings, test_freqs, test_phi, test_masses)
    mean_kinetic_eng_per_dof = mean(KE.(eachcol(velos), Ref(extended_test_masses))) / (3*N_atoms_test - 3)
    @test ustrip(u"eV", (0.5*test_kB*test_temp)) ≈ ustrip(u"eV",mean_kinetic_eng_per_dof) atol=1e-3

    _, velos = canonical_configs_and_velocities(settings, test_freqs, test_phi, test_masses)
    mean_kinetic_eng_per_dof = mean(KE.(eachcol(velos), Ref(extended_test_masses))) / (3*N_atoms_test - 3)
    @test ustrip(u"eV", (0.5*test_kB*test_temp)) ≈ ustrip(u"eV",mean_kinetic_eng_per_dof) atol=1e-3
end

@testset "QuantumHarmonicHeatCapacity" begin
    
end

