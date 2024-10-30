
function bose_einstein(freq, temp, kB, hbar)
    x =  upreferred((hbar * freq) / (kB * temp))
    return 1 / (exp(x) - 1)
end

function mean_amplitude(qc::QuantumConfigs, freq, mass)
    nᵢ = bose_einstein(freq, qc.temp, qc.kB, qc.hbar)
    return sqrt((qc.hbar * (2*nᵢ + 1)) / (2 * mass * freq))
end

function mean_amplitude(cc::ClassicalConfigs, freq, mass, temp)
    return sqrt((cc.kB*cc.temp)/mass) / freq
end

function check_mode(mode::Symbol)
    if mode ∉ [:quantum, :classical]
        throw(ArgumentError("Mode must be :quantum or :classical"))
    end
end

function extend_masses!(atom_masses, D)
    return collect(Iterators.flatten(Iterators.repeated(el, D) for el in atom_masses))
end

function prune_freqs_phi!(freqs, phi, D::Int)
    idx_rt = sortperm(abs.(freqs))
    if idx_rt != 1:D
        throw(ArgumentError("First D modes should be rigid translation modes, got $(idx_rt). There might be imaginary modes present"))
    end

    # Remove first D elements from freqs and phi
    freqs = freqs[D+1:end]
    phi = phi[:, D+1:end]

    return freqs, phi
end

#* Change signature to not enforce SVector/MVector
#* in examples use those though
function canonical_configs(CM::ConfigMode, freqs::AbstractVector{T},
                         phi::AbstractMatrix{T}, eq_positions::AbstractVector{SVector{D, Unitful.Length{T}}},
                         atom_masses::AbstractVector{Unitful.Mass{T}}, mode::Symbol;
                         nthreads::Int = Threads.nthreads()) where {T,D}

    check_mode(mode)

    N_atoms = length(eq_positions)

    # Remvoe rigid translation modes from freqs and phi
    # This makes it possible to broadcast without DivideByZero
    freqs, phi = prepare_freqs_phi!(freqs, phi, D)

    # Extend mass vector so there are D copies of each mass per atom
    # This makes broadcasting easier, could also rehspae phi to be (D*N_atoms - D) x N_atoms x 3
    extend_masses!(atom_masses, D)
    atom_masses_T = transpose(atom_masses) #& can this be done in place?

    # Create storage
    randn_storage = zeros(T, D*N_atoms)
    configs = zeros(eltype(first(eq_positions)), D*N_atoms, n_configs)

    time_unit = u"ps"
    velos = #*TODO

    # Pre-calculate phi * mean amplitudes
    phi_A = phi * mean_amplitude.(Ref(CM), freqs, atom_masses_T) # D*N_atoms - D x D*N_atoms 

    # bar = ProgressBar(1:n_configs; printing_delay = 0.1)
    # set_description(bar, "Making Configs")
    # p = Progress(n_configs; desc="Generating Configs", dt = 0.2)
    @tasks for n in 1:n_configs
        @set ntasks = nthreads
        randn!(randn_storage)
        configs[:, n] .= phi_A * randn_storage
        configs[:, n] .+= eq_positions
        # next!(p)
    end
    # finish!(p)

    return configs
end


# """
#     canonical_config!(CM::ConfigMode,
#                       config_storage::AbstractArray{Unitful.Length{T}, 1},
#                       randn_storage::AbstractVector{T},
#                       phi::AbstractMatrix{T},
#                       freqs::AbstractVector{Unitful.Frequency{T}},
#                       atom_masses_T::AbstractVector{Unitful.Mass{T}}
#                     ) where {T}

# Calculates a single canonical configuration.

# Parameters:
# -----------
# - `CM : ConfigMode`
#     Configuration data. Either `QuantumConfigSettings` or `ClassicalConfigSettings`
# - `config_storage : AbstractArray{Unitful.Length{T}, 1}`
#     Storage for the configuration. Typically a view to a larger array. Size : D*N_atoms x 1
# - `randn_storage : AbstractVector{T}`
#     Storage for the random numbers. Size : D*N_atoms x 1
# - `phi : AbstractMatrix{T}`
#     Eigenvectors of the dynamical matrix with the columns corresponding to the rigid
#     translation modes (ω = 0) removed. Size : D*N_atoms x D*N_atoms - D
# - `freqs : AbstractVector{Unitful.Frequency{T}}`
#     Frequencies of the normal modes with the rigid translation modes (ω = 0) removed.
#     Size : D*N_atoms - D x 1
# - `atom_masses_T : AbstractVector{Unitful.Mass{T}}`
#     Masses of the atoms as a row vector. Expects each mass to be duplicated `D` times per
#     atom to enable broadcasting. Size : 1 x D*N_atoms
# """
# function canonical_config!(CM::ConfigSettings,
#                             config_storage::AbstractArray{Unitful.Length{T}, 1},
#                             randn_storage::AbstractVector{T},
#                             phi::AbstractMatrix{T},
#                             freqs::AbstractVector{Unitful.Frequency{T}},
#                             atom_masses_T::AbstractVector{Unitful.Mass{T}}
#                         ) where {T}

#     randn!(randn_storage) #* spawn rng in parent task?

#     # Can I avoid the matrix allocation on amplitude?
#     config_storage .= phi_A * randn_storage

#     return config_storage
# end