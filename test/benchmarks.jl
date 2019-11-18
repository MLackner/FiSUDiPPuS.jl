using BenchmarkTools
using FiSUDiPPuS

options = joinpath(@__DIR__, "../data/default.jl")
include(options)

x,y = FiSUDiPPuS.get_data(settings) # settings is a dict in the options file

## Convert to correct types
start = settings[:datatype].(settings[:start])
bleach_weight = settings[:datatype](settings[:bleach_weight])

@btime FiSUDiPPuS.model($x[:,:], $start, $settings[:N], $bleach_weight);

x = 3000.0
A = [1.0, -1.0]
Γ = [8.0, 8.0]
ω = [2900.0, 2910.0]
@btime FiSUDiPPuS.sfspec($x, $A, $ω, $Γ)

# 496.331 μs (4505 allocations: 510.06 KiB)
# 336.573 μs (2259 allocations: 194.02 KiB)
# 296.531 μs (22 allocations: 89.84 KiB)
# 298.902 μs (19 allocations: 89.52 KiB)
# 196.338 μs (20 allocations: 89.63 KiB)


# 573.098 ns (11 allocations: 1.19 KiB)
# 447.626 ns (9 allocations: 928 bytes)
