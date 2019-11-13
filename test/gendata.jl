# Generate Data
A_ssp = [-4.0, 3.0, -8.0]
A_ppp = [8.0, -6.0, 3]
Γ     = [9.0,  8.0, 7.0]
ω     = [2870, 2910, 2940]
a_ssp = [1.0, 1.0, 1.0,
        0.99, 0.98, 0.97,
        0.98, 0.96, 0.94]
a_ppp = [1.0, 1.0, 1.0,
        0.96, 0.92, 0.88,
        0.98, 0.94, 0.90]


a_ssp = reshape(a_ssp, (3,3))' # num time steps, num resonances
a_ppp = reshape(a_ppp, (3,3))'
wn    = range(2830, stop=2990, length=301) |> collect

sig_ssp = zeros(size(a_ssp,1), length(wn))
sig_ppp = similar(sig_ssp)

for i = 1:size(a_ssp,1), j = 1:length(wn)
    sig_ssp[i,j] = FiSUDiPPuS.sfspec(wn[j], a_ssp[i,:] .* A_ssp, ω, Γ)
    sig_ppp[i,j] = FiSUDiPPuS.sfspec(wn[j], a_ppp[i,:] .* A_ppp, ω, Γ)
end

# Save Data
savepath = joinpath(@__DIR__, "../data")
d_ssp = Dict(
    "wavenumber" => wn,
    "signal"     => sig_ssp,
)
d_ppp = Dict(
    "wavenumber" => wn,
    "signal"     => sig_ppp,
)

save(joinpath(savepath, "sample_ssp.jld2"), d_ssp)
save(joinpath(savepath, "sample_ppp.jld2"), d_ppp)
