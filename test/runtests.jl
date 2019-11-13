using FiSUDiPPuS
using FileIO, JLD2
using Test

include("gendata.jl")

@testset "FiSUDiPPuS.jl" begin
    r = runfit(joinpath(@__DIR__, "../data/default.jl"))
    @test round(r.minimizer[1], digits=6) ≈ A_ssp[1]
    @test round(r.minimizer[4], digits=6) ≈ A_ppp[1]
end
