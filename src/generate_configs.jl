
function bose_einstein(freq, temp, kB, hbar)
    x =  upreferred((hbar * freq) / (kB * temp))
    return 1 / (exp(x) - 1)
end

function amplitude(qc::QuantumConfigs, freq, mass, temp)
    nᵢ = bose_einstein(freq, temp, qc.kB, qc.hbar)
    return sqrt((qc.hbar * (2*nᵢ + 1)) / (2 * mass * freq))
end

function amplitude(cc::ClassicalConfigs, freq, mass, temp)
    return sqrt((cc.kB*temp)/mass) / freq
end

function imaginary_mode_present(freqs::AbstractVector, tol = 1e-6)
    for freq in freqs
        if ustrip(abs(imag(freq))) > tol
            return true
        end
    end
    return false
end

# Assumes there are D modes with freuqency 0
function rigid_translation_modes(freqs, D::Int)
    idx_rt = sortperm(abs.(freqs))
    return SVector(idx_rt[1:D]...)
end

struct SelfConsistentConfigs{C,K,H,T}
    configs::Matrix{C}
    freq_checkpoints::Matrix{Float32}
    dynmat_checkpoints::Array{Float32, 3}
    checkpoint_idxs::Vector{Int}
    kB::K
    hbar::H
    temp::T
    n_iters::Int
end

function check_mode(mode::Symbol)
    if mode ∉ [:quantum, :classical]
        throw(ArgumentError("Mode must be :quantum or :classical"))
    end
end

function extend_masses!(atom_masses, D)

end

function prune_rt_modes!(freqs, phi, rtm_idxs)
    # Remove rigid translation modes from freqs and phi
    freqs = deleteat!(freqs, rtm_idxs)
end

#* Change signature to not enforce SVector/MVector
#* in examples use those though
function canonical_configs(CM::ConfigMode, N_atoms::Int, freqs::AbstractVector{T},
                         phi::AbstractMatrix{T}, eq_positions::AbstractVector{SVector{D, Unitful.Length{T}}},
                         atom_masses::AbstractVector{Unitful.Mass{T}}, temp::Unitful.Temperature{T}, mode::Symbol;
                         nthreads::Int = Threads.nthreads()) where {T,D}

    check_mode(mode)

    # Remvoe rigid translation modes from freqs and phi
    rtm_idxs = rigid_translation_modes(freqs, D)
    prune_rt_modes!(freqs, phi, rtm_idxs)

    # Extend mass vector so there are D copies of each mass per atom
    # This makes broadcasting over N_dof easier
    extend_masses!(atom_masses, D)
    atom_masses_T = transpose(atom_masses) # can this be done in place?

    # Create storage
    randn_storage = zeros(T, N_dof - D)
    # config_storage = zeros(MVector{D, Unitful.Length{T}}, N_dof, n_configs) #* idk how to do this with mvectors

    # bar = ProgressBar(1:n_configs; printing_delay = 0.1)
    # set_description(bar, "Making Configs")
    # p = Progress(n_configs; desc="Generating Configs", dt = 0.2)
    @tasks for n in 1:n_configs
        @set ntasks = nthreads
        canonical_config!(config_storage[n], randn_storage, phi, freqs, atom_masses_T)
        configs[:, n] .+= eq_positions
        # next!(p)
    end
    # finish!(p)

    return configs
end


function canonical_config!(
                            config_storage::AbstractArray{MVector{D, Unitful.Length{T}}, 1}, #typically passed as view to larger array
                            randn_storage::AbstractMatrix{T}, #unitless
                            phi::AbstractMatrix{T}, #unitless
                            freqs::AbstractVector{Unitful.Frequency{T}},
                            atom_masses_T::AbstractVector{Unitful.Mass{T}}
                        ) where {D,T}
    # Assumes `freqs` are passed with the rigid translation modes removed
        # to facilitate broadcasting. The corresponding columns from `phi`` should also be removed.
    
    # Assumes mass vector is duplicated  so that each dof has a mass
        # Again to facilitate broadcasting. Also this should be Transposed.

    randn!(randn_storage) #* spawn rng in parent task?

    # amplitude(freqs[m], atom_masses[i], kB, temp) --> should be N_dof x len(atom_masses)
    # configs[ii] += (A * z[m] * phi[ii, m])

    #& think I need to broadcast over masses transposed
    config_storage .= amplitude.(CM, freqs, atom_masses_T, temp) .* (phi * randn_storage) #* dont think this works? mass/freq indexing weird

    return config_storage
end