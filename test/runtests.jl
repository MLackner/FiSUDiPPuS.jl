using FiSUDiPPuS
using FileIO, JLD2
using BlackBoxOptim
using Test
using FiSUDiPPuS

@testset "Parameter Decomposing" begin
    A_in = 4.0
    χnr_in = 1.5
    Δω_in = 2.0
    β_in  = 3.0
    ω_in = 2900.0
    Γ_in = 10.0
    p = [A_in, χnr_in, Δω_in, β_in, ω_in, Γ_in]

    A, φ, ω, Γ, a, δω, χnr, Δω, β = FiSUDiPPuS.decompose_parameters(p)

    @test A == [A_in]
    @test φ == [0.0]
    @test ω == [ω_in]
    @test Γ == [Γ_in]
    @test a == [1.0]
    @test δω == [0.0]
    @test χnr == χnr_in
    @test Δω == Δω_in
    @test β == 1.0 # for tstep=0 no β should stay neutral

    A_in = [1.0, 2.0]
    φ_in = [π, -π]
    χnr_in = 2.0
    Δω_in = 3.0
    β_in = 4.0
    ω_in = [2890.0, 2910.0]
    Γ_in = [5.0, 6.0]
    a_in = [0.99, 0.99, 0.9, 0.8]
    δω_in = [1.0, 1.0, -0.1, -0.2]

    p = [
        A_in...,
        φ_in...,
        χnr_in...,
        Δω_in...,
        β_in...,
        ω_in...,
        Γ_in...,
        a_in...,
        δω_in...,
    ]
    A, φ, ω, Γ, a, δω, χnr, Δω, β = FiSUDiPPuS.decompose_parameters(p;
            N=2, phase=true, ω_shift=true, n_steps=2, tstep=1
        )

    @test A == A_in
    @test φ == φ_in
    @test ω == ω_in
    @test Γ == Γ_in
    @test a == a_in[1:2]
    @test δω == δω_in[1:2]
    @test χnr == χnr_in
    @test Δω == Δω_in
    @test β == β_in

    A, φ, ω, Γ, a, δω, χnr, Δω, β = FiSUDiPPuS.decompose_parameters(p;
            N=2, phase=true, ω_shift=true, n_steps=2, tstep=2
        )
    @test a == a_in[3:4]
    @test δω == δω_in[3:4]

    A_in = [1.0, 2.0]
    φ_in = [π, -π]
    χnr_in = 2.0
    Δω_in = 3.0
    β_in = 4.0
    ω_in = [2890.0, 2910.0]
    Γ_in = [5.0, 6.0]
    a_in = [1.0, 1.0, 0.9, 0.8]

    p = [
        A_in...,
        φ_in...,
        χnr_in...,
        Δω_in...,
        β_in...,
        ω_in...,
        Γ_in...,
        a_in...,
    ]
    A, φ, ω, Γ, a, δω, χnr, Δω, β = FiSUDiPPuS.decompose_parameters(p;
            N=2, phase=true, ω_shift=false, n_steps=2, tstep=1
        )

    @test A == A_in
    @test φ == φ_in
    @test ω == ω_in
    @test Γ == Γ_in
    @test a == a_in[1:2]
    @test δω == [0.0, 0.0]
    @test χnr == χnr_in
    @test Δω == Δω_in
    @test β == β_in

    A_in = 4.0
    χnr_in = 1.5
    Δω_in = 2.0
    β_in  = 3.0
    ω_in = 2900.0
    Γ_in = 10.0
    a_in = 0.5
    p = [A_in, χnr_in, Δω_in, β_in, ω_in, Γ_in, a_in, ]

    @show A, φ, ω, Γ, a, δω, χnr, Δω, β = FiSUDiPPuS.decompose_parameters(p;
        phase=false, tstep=1, n_steps=1, N=1, ω_shift=false,
    )

    @test A == [A_in]
    @test φ == [0.0]
    @test ω == [ω_in]
    @test Γ == [Γ_in]
    @test a == [a_in]
    @test δω == [0.0]
    @test χnr == χnr_in
    @test Δω == Δω_in
    @test β == β_in
end

@testset "Model" begin
    A = 1.0
    φ = 0.0
    χnr = 0.0
    Δω = 0.0
    β = 1.0
    ω = 2900
    Γ = 1.0

    x = 2800:3000

    p = [
        A,
        χnr,
        Δω,
        β,
        ω,
        Γ
    ]

    y = model(x,p; phase=false, tstep=0, n_steps=0, N=1, ω_shift=false)

    @test maximum(y) == 1.0
    @test x[argmax(y)] == 2900

    p[3] = -2.0 # change Δω

    y = model(x,p; phase=false, tstep=0, n_steps=0, N=1, ω_shift=false)

    @test maximum(y) == 1.0
    @test x[argmax(y)] == 2898

    a = 0.5

    p = [
        A,
        χnr,
        Δω,
        β,
        ω,
        Γ,
        a,
    ]

    y = model(x,p; phase=false, tstep=0, n_steps=1, N=1, ω_shift=false)

    @test maximum(y) == 1.0
    @test x[argmax(y)] == 2900

    y = model(x,p; phase=false, tstep=1, n_steps=1, N=1, ω_shift=false)

    @test minimum(y) == -0.75
    @test x[argmin(y)] == 2900

end

## MAKE A SPECTRUM

A1 = [0.5, 2.0, -0.3]
A2 = [-1.0, 0.0, 0.3]
a = [0.5, 0.9, 0.9]
χnr = 0.0
Δω1 = 0.0
Δω2 = 1.5
β1 = 1.0
β2 = 1.5
ω = [2923, 2938, 2964]
Γ = [9.0, 9.0, 9.0]

p = [A1..., χnr, Δω1, β1,
     A2..., χnr, Δω2, β2,
     ω..., Γ..., a...]

x = 2800:3000
n_steps = 1
N = 3

signal1 = zeros(n_steps+1,length(x))
signal2 = zeros(n_steps+1,length(x))
for i in 0:n_steps
    signal1[i+1,:] = model(x,p;
        phase=false, n_steps=n_steps, N=N,
        ω_shift=false, tstep=i, specidx=1
    )
    signal2[i+1,:] = model(x,p;
        phase=false, n_steps=n_steps, N=N,
        ω_shift=false, tstep=i, specidx=2
    )
end

s1 = Spectrum(x, signal1)
s2 = Spectrum(x, signal2)

spectra = [s1, s2]
settings = Dict(
    :options => Dict(
        :SearchRange => Tuple{Float64,Float64}[
            # SPEC 1
            # oscillator strengths
            fill((-1.5,3), N)...,
            (0,0),  # χnr
            (0,0), # Δω
            (1,1), # β
            # SPEC 2
            # oscillator strengths
            fill((-1.5,3), N)...,
            (0,0),  # χnr
            (-2,2), # Δω
            (0.5,2.0), # β
            # UNIVERSAL
            # wavenumbers
            (2910, 2930),
            (2930, 2950),
            (2950, 2970),
            # damping coefficients
            fill((6.0, 15.0), N)...,
            # a values
            fill((0.3,1.1), N)...,
        ],
        :MaxSteps => 1e5,
        :MaxTime => 20,
    ),
    :n_steps => 1,
    :ω_shift => false,
    :phase => false,
    :N => 3,
    :bleach_weight => 3,
)
result = runfit(spectra, settings)

# plt = plot()
# for i in 1:size(signal,1)
#     plot!(x, signal[i,:])
# end
# display(plt)
