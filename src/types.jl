abstract type ConfigSettings end

@derived_dimension BoltzmannConstUnits Unitful.𝐌*(Unitful.𝐋^2)*(Unitful.𝐓^-2)*(Unitful.𝚯^-1) true
@derived_dimension MolarBoltzmannConstUnits Unitful.𝐌*(Unitful.𝐋^2)*(Unitful.𝐓^-2)*(Unitful.𝚯^-1)*(Unitful.𝐍^-1) true

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