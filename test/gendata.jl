# Generate Data
A_ssp = [-4.0, 3.0, -8.0]
A_ppp = [8.0, -6.0, 3]
Γ     = [9.0,  8.0, 7.0]
ω     = [2870, 2910, 2940]
a_ssp = [1.0, 1.0, 1.0,
        0.91, 0.92, 0.93,
        0.90, 0.90, 0.90]
a_pow = 0.9   # factor by which we power the a_ssp values to
               # to get the a_ppp values
χ3 = 0.1
φ  = -π

a_ssp = reshape(a_ssp, (3,3))' # num time steps, num resonances
a_ppp = a_ssp .^ a_pow
wn    = range(2830, stop=2990, length=301) |> collect

sig_ssp = zeros(size(a_ssp,1), length(wn))
sig_ppp = similar(sig_ssp)

for i = 1:size(a_ssp,1), j = 1:length(wn)
    sig_ssp[i,j] = FiSUDiPPuS.sfspec(wn[j], a_ssp[i,:] .* A_ssp, ω, Γ, χ3=χ3, φ=φ)
    sig_ppp[i,j] = FiSUDiPPuS.sfspec(wn[j], a_ppp[i,:] .* A_ppp, ω, Γ, χ3=χ3, φ=φ)
end

# Save Data
savepath = joinpath(@__DIR__, "../data")
d_ssp = Dict(
    "wavenumber" => wn,
    "signal"     => sig_ssp,
    "dltime"     => range(-3, 5, length=length(a_ssp)÷3)
)
d_ppp = Dict(
    "wavenumber" => wn,
    "signal"     => sig_ppp,
    "dltime"     => range(-3, 5, length=length(a_ssp)÷3)
)

save(joinpath(savepath, "sample_ssp.jld2"), d_ssp)
save(joinpath(savepath, "sample_ppp.jld2"), d_ppp)
