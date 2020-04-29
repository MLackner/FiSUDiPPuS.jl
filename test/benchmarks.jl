using BenchmarkTools
using FiSUDiPPuS

N = 3
n_spec = 1
n_steps = 2

# A parameters and χ, Δω, β
A = rand(N*n_spec + 3*n_spec)
ω = rand(2800:3000, N)
Γ = rand(N)
a = rand(N * n_steps + N)

p = [A..., ω..., Γ..., a...]

x = range(2800, 3000, length=512)

A, ω, Γ, φ, χnr = (rand(5), rand(5), rand(5), fill(0.0, 5), 0.0)
Γc = Γ .|> ComplexF64
φc = φ .|> ComplexF64
Ac = A .|> ComplexF64
ωc = ω .|> ComplexF64
χnrc = χnr |> ComplexF64
@btime FiSUDiPPuS.sfspec(x[1], A, ω, Γ, φ; χnr = χnr)
# 537.707 ns (4 allocations: 64 bytes)
@btime FiSUDiPPuS.sfspec(x[1], Ac, ωc, Γc, φc; χnr = χnrc)
# 509.188 ns (4 allocations: 64 bytes)
# 486.272 ns (4 allocations: 64 bytes)
@btime FiSUDiPPuS.sfspec(x[1], Ac, ωc, Γc)
# 116.253 ns (2 allocations: 32 bytes)
# I can't get this four times speedup translated to the model function
@btime FiSUDiPPuS.decompose_parameters(p;
    n_steps=n_steps,
    N=N,
    specidx=1,
    phase=false,
    ω_shift=false,
    tstep=0
)
# 967.600 ns (14 allocations: 1.45 KiB)
@btime FiSUDiPPuS.model(x,p;
    phase=false,
    ω_shift=false,
    n_steps=n_steps,
    N=N,
    tstep=0
)
# 91.503 μs (1040 allocations: 117.63 KiB)
# 58.314 μs (1056 allocations: 119.17 KiB)
# 30.020 μs (32 allocations: 7.14 KiB)
# 24.292 μs (38 allocations: 7.91 KiB)
# 20.660 μs (38 allocations: 7.88 KiB)
@btime FiSUDiPPuS.model(x,p;
    phase=false,
    ω_shift=false,
    n_steps=n_steps,
    N=N,
    tstep=1
)
# 186.784 μs (1553 allocations: 173.73 KiB)
# 136.760 μs (1569 allocations: 175.28 KiB)
# 117.164 μs (1584 allocations: 176.80 KiB)
# 57.173 μs (48 allocations: 8.75 KiB)
# 45.478 μs (54 allocations: 9.53 KiB)
# 38.661 μs (54 allocations: 9.47 KiB)
