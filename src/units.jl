# When energy is for exmaple, kcal/mol, the freuency squared has units kcal / mol / Å^2
# @derived_dimension MolarFrequency (𝐍^-1/2) * inv(𝐓) true

@derived_dimension BoltzmannConstUnits Unitful.𝐌*(Unitful.𝐋^2)*(Unitful.𝐓^-2)*(Unitful.𝚯^-1) true
@derived_dimension MolarBoltzmannConstUnits Unitful.𝐌*(Unitful.𝐋^2)*(Unitful.𝐓^-2)*(Unitful.𝚯^-1)*(Unitful.𝐍^-1) true

@derived_dimension hBarUnits (Unitful.𝐋^2)* Unitful.𝐌 * (Unitful.𝐓^-1) true
@derived_dimension MolarhBarUnits (Unitful.𝐋^2)* Unitful.𝐌 * (Unitful.𝐓^-1) * (Unitful.𝐍^-1) true

function check_units(kB, hbar)
    both_molar = (kB isa MolarBoltzmannConstUnits) && (hbar isa MolarhBarUnits)
    both_not_molar = (kB isa BoltzmannConstUnits) && (hbar isa hBarUnits)

    if both_molar || both_not_molar
        return true
    else
        return false
    end
end