export QuantumConfigSettings, ClassicalConfigSettings

abstract type ConfigSettings end

const DefaultFloat = Float32

struct QuantumConfigSettings{K,H} <: ConfigSettings 
    kB::K
    h_bar::H
    n_configs::Int
    temperature::typeof(DefaultFloat(1.0u"K"))
end

function QuantumConfigSettings(kB, h_bar, n_configs, temperature)
    if check_units(kB, h_bar)
        return QuantumConfigSettings{typeof(kB), typeof(h_bar)}(
                        DefaultFloat(kB), DefaultFloat(h_bar), n_configs, DefaultFloat(temperature))
    else
        throw(ArgumentError("Units of kB and h_bar are not comensurate. Must both be molar or non-molar."))
    end
end

struct ClassicalConfigSettings{K} <: ConfigSettings
    kB::K
    n_configs::Int
    temperature::typeof(DefaultFloat(1.0u"K"))
end

function ClassicalConfigSettings(kB, n_configs, temperature)
    if kB isa BoltzmannConstUnits || kB isa MolarBoltzmannConstUnits
        return ClassicalConfigSettings{typeof(kB)}(DefaultFloat(kB), n_configs, DefaultFloat(temperature))
    else
        throw(ArgumentError("kB does not have proper units."))
    end
end