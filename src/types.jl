abstract type ConfigSettings end

@derived_dimension BoltzmannConstUnits 𝐌*(𝐋^2)*(𝐓^-2)*(𝚯^-1) true
@derived_dimension MolarBoltzmannConstUnits 𝐌*(𝐋^2)*(𝐓^-2)*(𝚯^-1)*(𝐍^-1) true

@derived_dimension hBarUnits (𝐋^2)* 𝐌 * (𝐓^-1)
@derived_dimension MolarhBarUnits

struct QuantumConfigSettings{T} <: ConfigSettings 
    kB::Union{BoltzmannConstUnits{T}, MolarBoltzmannConstUnits{T}}
    h_bar::Union{hBarUnits{T}, MolarhBarUnits{T}}
    n_configs::Int
    temperature::Unitful.Tempearture{T}
end

struct ClassicalConfigSettings{T} <: ConfigSettings
    kB::Union{BoltzmannConstUnits{T}, MolarBoltzmannConstUnits{T}}
    n_configs::Int
    tempearture::Unitful.Temperature{T}
end