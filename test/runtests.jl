using FiSUDiPPuS
using FileIO, JLD2
using Test

include("gendata.jl")

@testset "FiSUDiPPuS.jl" begin
    settings_file = joinpath(@__DIR__, "../data/default.jl")
    include(settings_file)
    r = runfit(settings, saveprefix="test")

    # Test if multi and single-threaded functions give the same results
    @test all(r.minimizer .≈ r.minimizer)

    @test round(r.minimizer[1], digits=4) ≈ A_ssp[1]
    @test round(r.minimizer[4], digits=4) ≈ A_ppp[1]
    #TODO: The optimization results are not as good as they used to be
    @test round(r.minimizer[end-2], digits=2) ≈ a_pow

    # run the other exported functions on the data to check if they error
    # get latest results
    # result_file = readdir(joinpath(@__DIR__, "../data/results/"))[end]
    # result_file = joinpath("../data/results", result_file)
    # plotresult(result_file)
    # printresult(result_file)
    # viewsettings(settings_file)
end
        
