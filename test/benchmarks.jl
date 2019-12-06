using BenchmarkTools
using FiSUDiPPuS

p = [
    1.0, -2.0,  # A_ssp
    2.0, -0.5,  # A_ppp
    2870, 2940, # ω
    10.0, 10.0, # Γ
    0.9, 0.95,  # a
    1.0,        # Δω_ppp
    1.5,        # a_pow
    0.02,        # χ3
    -π,          # φ
]

x = range(2800, 3050, length=151) |> collect
# y = model(x,p, diff=false, pol=:ppp)
# plot(x,y)

# Preallocate complex array
y = zeros(length(x))
@btime model!(y, x, p, diff=false, pol=:ssp)
@btime model!(y, x, p, diff=true , pol=:ssp)
@btime model!(y, x, p, diff=false, pol=:ppp)
@btime model!(y, x, p, diff=true , pol=:ppp)

# 11.742 μs (8 allocations: 5.56 KiB)
# 23.296 μs (12 allocations: 9.58 KiB)
# 11.786 μs (9 allocations: 5.66 KiB)
# 23.498 μs (14 allocations: 9.77 KiB)

# 11.420 μs (6 allocations: 1.73 KiB)
# 22.712 μs (9 allocations: 3.25 KiB)
# 11.450 μs (7 allocations: 1.83 KiB)
# 22.934 μs (11 allocations: 3.44 KiB)

# 10.790 μs (5 allocations: 224 bytes)
# 28.772 μs (157 allocations: 14.42 KiB)
# 16.843 μs (156 allocations: 14.38 KiB)
# 48.191 μs (459 allocations: 42.73 KiB)

# 10.102 μs (6 allocations: 512 bytes)
# 20.101 μs (7 allocations: 608 bytes)
# 10.100 μs (6 allocations: 512 bytes)
# 20.093 μs (7 allocations: 608 bytes)
