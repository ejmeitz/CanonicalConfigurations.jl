# CanonicalConfigurations.jl
 Generates samlpes from the classical or quantum canonical ensemble from phonon modes.


# Example:

```julia
using Unitful
using CanonicalConfigurations

temp = 300u"K"
kB = uconvert(u"kcal / mol / K", Unitful.k*Unitful.Na)

masses = 39.95u"g/mol" * ones(N_atoms)

## Load frequencies and eigenvectors
## On you to calculate them


freq_unit = sqrt(u"kcal / Å^2 / g")
freqs *= freq_unit

n_configs = 250_000
settings = ClassicalConfigSettings(kB, n_configs, temp)
configs = canonical_configs(settings, freqs, phi, masses)
```
