using Random
using FileIO, JLD2
Random.seed!(1)

# Generate Data
A_ssp = [-2.0, -4.0, 8.0]
A_ppp = [-9.0, 4.0, 1.0]
Γ     = [10.0, 8.0, 9.0]
ω     = [2880, 2910, 2930]
a_ssp = [1.0, 1.0, 1.0,
        0.85, 0.90, 0.91]
δω    = [0.0, 0.0, 0.0,
         0.0, -2.0, 0.0]
Δω = -0.0
a_pow = 0.9   # factor by which we power the a_ssp values to
               # to get the a_ppp values
χ3 = 0.1
φ  = π/2

a_ssp = reshape(a_ssp, (3,2))' #  num resonances, num time steps
δω    = reshape(δω   , (3,2))' #  num resonances, num time steps
a_ppp = a_ssp .^ a_pow
wn    = range(2830, stop=2990, length=301) |> collect

sig_ssp = zeros(size(a_ssp,1), length(wn))
sig_ppp = similar(sig_ssp)

for i = 1:size(a_ssp,1), j = 1:length(wn)
    sig_ssp[i,j] = FiSUDiPPuS.sfspec(wn[j], a_ssp[i,:] .* A_ssp, ω .+ δω[i,:], Γ, χ3=χ3, φ=φ)
    sig_ppp[i,j] = FiSUDiPPuS.sfspec(wn[j], a_ppp[i,:] .* A_ppp, ω .+ Δω .+ δω[i,:], Γ, χ3=χ3, φ=φ)
end

# sig_ssp .+= randn(size(sig_ssp)) ./ 200
# sig_ppp .+= randn(size(sig_ppp)) ./ 200

# Save Data
savepath = joinpath(@__DIR__, "../data")
d_ssp = Dict(
    "wavenumber" => wn,
    "signal"     => sig_ssp,
    "dltime"     => range(-3, 5, length=size(a_ssp,1))
)
d_ppp = Dict(
    "wavenumber" => wn,
    "signal"     => sig_ppp,
    "dltime"     => range(-3, 5, length=size(a_ssp,1))
)

save(joinpath(savepath, "sample_ssp.jld2"), d_ssp)
save(joinpath(savepath, "sample_ppp.jld2"), d_ppp)
