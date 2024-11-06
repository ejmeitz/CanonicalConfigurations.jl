using Pkg
Pkg.activate("C:/Users/ejmei/.julia/environments/research")

using SimpleCrystals 
using JLD2
using Unitful
using ForceConstants
using CanonicalConfigurations
using Statistics
using ProgressMeter
using OhMyThreads
using Plots

function KE(v, m)
    return 0.5 * sum(v.^2 .* m)
end

function energy_from_potential(all_configs, sys_eq::SuperCellSystem, pot::Potential, energy_unit)
    n_configs = size(all_configs,2)
    N_atoms = n_atoms(sys_eq)

    U = zeros(Float32, n_configs) * energy_unit
    L_unit = unit(first(all_configs))

    p = Progress(size(all_configs,2); desc="Calculating Potential", dt = 0.2, color = :blue)
    @tasks for j in 1:n_configs
        @set begin
            ntasks = Threads.nthreads()
            scheduler = :static
        end
        @local posns = [MVector{3}(0.0, 0.0, 0.0)*L_unit for _ in 1:N_atoms]
        @views config = all_configs[:,j]

        #Update posn vectors
        for i in 1:N_atoms
            # Add Equilibrium positions
            @views posns[i] .= config[3*(i-1) + 1 : 3*(i-1) + 3] .+ sys_eq.atoms[i].position
        end

        # posns = [SVector{3}(config[3*(i-1) + 1], config[3*(i-1) + 2], config[3*(i-1) + 3]) for i in 1:N_atoms]
        U[j] = energy_loop(pot, posns, sys_eq.box_sizes_SC, N_atoms, pot.r_cut)
        next!(p)
    end
    finish!(p)
    
    return U
end

function run_loop_var(ifc_path, pot::StillingerWeberSilicon, crys::Crystal, temp, n_configs::Int)
    @info "Running for $(temp)K"
    freqs_sq, phi = load(ifc_path, "freqs_sq", "phi")
    freqs = Float32.(sqrt.(freqs_sq))
    phi = Float32.(phi)

    freq_unit = sqrt(u"eV / Å^2 / u") #* ASSUME SW
    freqs *= freq_unit
    atom_masses = Float32.(atomic_mass(crys))

    
    kB = uconvert(u"eV / K", Unitful.k)
    h_bar = uconvert(u"eV*s", Unitful.ħ)
    settings = QuantumConfigSettings(kB, h_bar, n_configs, temp)
    # settings = ClassicalConfigSettings(kB, n_configs, temp)
    configs, velos = canonical_configs_and_velocities(settings, freqs, phi, atom_masses)

    U = energy_from_potential(configs, sys_eq, pot, energy_unit(pot))
    extended_masses = collect(Iterators.flatten(Iterators.repeated(el, 3) for el in atom_masses))
    K = KE.(eachcol(velos), Ref(extended_masses))
    E = U .+ K
    cv = upreferred(var(E) / (temp * temp * kB) / (kB * 3* n_atoms(sys_eq)))
    @show cv
    cv_u = upreferred(var(U) / (temp * temp * kB) / (kB * 3* n_atoms(sys_eq)))
    cv_K = upreferred(var(K) / (temp * temp * kB) / (kB * 3* n_atoms(sys_eq)))
    @show cv_u
    @show cv_K

    return cv
end

function regression_approx(avg_energies, avg_energy_std, temps)
    # x_observed = measurement.(temps, 0.0)
    # y_observed = measurement.(avg_energies, avg_energy_std)
    linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
    (i, s) = linreg(temps, avg_energies)
    return [i,s] #allocates
end

function run_loop_dUdT(ifc_path, pot::StillingerWeberSilicon, crys::Crystal, base_temp, n_configs::Int, dTs_percents)
    @info "Running for $(base_temp)K"
    freqs_sq, phi = load(ifc_path, "freqs_sq", "phi")
    freqs = Float32.(sqrt.(freqs_sq))
    phi = Float32.(phi)

    freq_unit = sqrt(u"eV / Å^2 / u") #* ASSUME SW
    freqs *= freq_unit
    atom_masses = Float32.(atomic_mass(crys))

    
    kB = uconvert(u"eV / K", Unitful.k)
    h_bar = uconvert(u"eV*s", Unitful.ħ)

    dTs = base_temp .* dTs_percents
    
    mean_Us = []
    mean_Ks = []
    mean_Es = []
    U_stds = []
    K_stds = []
    E_stds = []
    for (i,dT) in enumerate(dTs)
        settings = QuantumConfigSettings(kB, h_bar, n_configs, base_temp + dT)
        # settings = ClassicalConfigSettings(kB, n_configs, base_temp + dT)
        configs, velos = canonical_configs_and_velocities(settings, freqs, phi, atom_masses)

        U = energy_from_potential(configs, sys_eq, pot, energy_unit(pot))
        extended_masses = collect(Iterators.flatten(Iterators.repeated(el, 3) for el in atom_masses))
        
        K = KE.(eachcol(velos), Ref(extended_masses))
        E = U .+ K

        push!(mean_Us, mean(U))
        push!(mean_Ks, mean(K))
        push!(mean_Es, mean(E))
        push!(U_stds, std(U))
        push!(K_stds, std(K))
        push!(E_stds, std(E))
    end

    
    mean_Us = upreferred.(mean_Us)
    mean_Ks = upreferred.(mean_Ks)
    mean_Es = upreferred.(mean_Es)
    U_stds = upreferred.(U_stds)
    K_stds = upreferred.(K_stds)
    E_stds = upreferred.(E_stds)

    kB = upreferred(kB)

    temps = base_temp .+ dTs
    _, slope_U = regression_approx(ustrip.(mean_Us), ustrip.(U_stds), ustrip.(temps))
    slope_U *= unit(first(mean_Us)) / unit(first(temps))
    cv_U = slope_U ./ (3* n_atoms(sys_eq) * kB)

    _, slope_K = regression_approx(ustrip.(mean_Ks), ustrip.(K_stds), ustrip.(temps))
    slope_K *= unit(first(mean_Ks)) / unit(first(temps))
    cv_K = slope_K ./ (3* n_atoms(sys_eq) * kB)

    _, slope_E = regression_approx(ustrip.(mean_Es), ustrip.(E_stds), ustrip.(temps))
    slope_E *= unit(first(mean_Es)) / unit(first(temps))
    cv_E = slope_E ./ (3* n_atoms(sys_eq) * kB)

    @show cv_U
    @show cv_K
    @show cv_E

    return cv_U, cv_K, cv_E
end

function bose_einstein(freq, temp, kB, hbar)
    x =  upreferred((hbar * freq) / (kB * temp))
    return 1 / (exp(x) - 1)
end

function hld_heat_capapacity(kB, h_bar, T, freqs)
    nᵢ = bose_einstein.(freqs, T, kB, h_bar)
    return sum((1/(kB*T*T)) .* (h_bar * freqs).^2 .* nᵢ .* (nᵢ .+ 1))
end


crys = Diamond(5.43u"Å", :Si, SVector(3,3,3))
pot = StillingerWeberSilicon()
sys_eq = SuperCellSystem(crys)

n_configs = 40000
temps = [100, 300, 500, 700, 900, 1100, 1300] * u"K"
path = (T) -> "Z:/emeitz/Data/ForceConstants/AvgINM_SW/AvgIFC_SW_3UC_$(ustrip(T))K_CLEANED.jld2"

paths = path.(temps)

# cvs = run_loop_var.(paths, Ref(pot), Ref(crys), temps, Ref(n_configs))

dTs_percents = [-0.01, -0.005, 0.0, 0.005, 0.01]
heat_caps = run_loop_dUdT.(paths, Ref(pot), Ref(crys), temps, n_configs, Ref(dTs_percents))

hld_heat_caps = []
kB = uconvert(u"eV / K", Unitful.k)
h_bar = uconvert(u"eV*s", Unitful.ħ)
for temp in temps
    freqs_sq = load(path(temp), "freqs_sq")
    freqs = Float32.(sqrt.(freqs_sq))

    freqs = freqs[4:end] #remove rigid translation

    freq_unit = sqrt(u"eV / Å^2 / u") #* ASSUME SW
    freqs *= freq_unit

    push!(hld_heat_caps, hld_heat_capapacity(kB, h_bar, temp, freqs))
end

cv_Us = []
cv_Ks = []
cv_Es = []
for (cv_U, cv_K, cv_E) in heat_caps
    push!(cv_Us, cv_U)
    push!(cv_Ks, cv_K)
    push!(cv_Es, cv_E)
end

scatter(temps, cv_Es, label="d<E>/dT", xlabel="Temperature",
         ylabel="Cv / N kB", title="Heat Capacity of SW Silicon (Canonical Configs)",
         legend = :bottomright, markersize = 7)
scatter!(temps, upreferred.(hld_heat_caps /(3*n_atoms(sys_eq) * kB)), 
            label="Harmonic Lattice Dynamics", markersize = 7)
savefig("C:/Users/ejmei/Box/Research/Projects/QuantumCv/SW_Cv_vs_HLD.png")