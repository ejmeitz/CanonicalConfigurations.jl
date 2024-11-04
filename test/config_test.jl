using SimpleCrystals 
using JLD2
using Unitful

crys = Diamond(5.43u"Å", :Si, SVector(3,3,3))

kB = uconvert(u"eV / K", Unitful.k)
n_configs = 2000
temp = 1300u"K"
settings = ClassicalConfigSettings(kB, n_configs, temp)

path = "Z:/emeitz/Data/ForceConstants/AvgINM_SW/AvgIFC_SW_3UC_$(ustrip(temp))K_CLEANED.jld2"
freqs_sq, phi = load(path, "freqs_sq", "phi")
freqs = Float32.(sqrt.(freqs_sq))
phi = Float32.(phi)

freq_unit = sqrt(u"eV / Å^2 / u")
freqs *= freq_unit

# eq_positions = position(crys) #* make svector?
atom_masses = Float32.(atomic_mass(crys))

canonical_configs(settings, freqs, phi, atom_masses)
canonical_configs_and_velocities(settings, freqs, phi, atom_masses)