using FiSUDiPPuS
using FileIO, JLD2
using Test

include("gendata.jl")

@testset "FiSUDiPPuS.jl" begin
    settings_file = joinpath(@__DIR__, "../data/default.jl")
    r = runfit(settings_file; saveprefix="test")
    @test round(r.minimizer[1], digits=6) ≈ A_ssp[1]
    @test round(r.minimizer[4], digits=6) ≈ A_ppp[1]
    @test round(r.minimizer[end-1], digits=6) ≈ a_pow

    # run the other exported functions on the data to check if they error
    # get latest results
    result_file = readdir(joinpath(@__DIR__, "../data/results/"))[end]
    result_file = joinpath("../data/results", result_file)
    plotresult(result_file)
    printresult(result_file)
    viewsettings(settings_file)
end
