module CanonicalConfigurations


using Unitful
using OhMyThreads
using ThreadPinning
using ProgressMeter
using Random

pinthreads(:cores)

include("types.jl")
include("units.jl")
include("generate_configs.jl")
include("io.jl")

end