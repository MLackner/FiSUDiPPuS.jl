using FiSUDiPPuS
using FileIO, JLD2
using BlackBoxOptim
using Test

include("gendata.jl")

@testset "FiSUDiPPuS.jl" begin
    settings_file = joinpath(@__DIR__, "../data/default.jl")
    include(settings_file)
    r = runfit(settings, saveprefix="test")

    p = best_candidate(r)
    @test round(p[1], digits=1) * 10 ≈ A_ssp[1]
    @test round(p[4], digits=1) * 10 ≈ A_ppp[1]
    #TODO: The optimization results are not as good as they used to be
    @test round(p[end-2], digits=1) ≈ a_pow

    # run the other exported functions on the data to check if they error
    # get latest results
    # result_file = readdir(joinpath(@__DIR__, "../data/results/"))[end]
    # result_file = joinpath("../data/results", result_file)
    # plotresult(result_file)
    # printresult(result_file)
    # viewsettings(settings_file)
end

# begin
#     settings_file = joinpath(@__DIR__, "../data/default.jl")
#     include(settings_file)
#     algs = [
#         # Natural Evolution Strategies:
#         # Separable NES: separable_nes
#         # Exponential NES:
#         :xnes,
#         # Distance-weighted Exponential NES:
#         :dxnes,
#         # Differential Evolution optimizers, 5 different:
#         # Adaptive DE/rand/1/bin:
#         :adaptive_de_rand_1_bin,
#         # Adaptive DE/rand/1/bin with radius limited sampling:
#         :adaptive_de_rand_1_bin_radiuslimited,
#         # DE/rand/1/bin:
#         :de_rand_1_bin,
#         # DE/rand/1/bin with radius limited sampling (a type of trivial geography):
#         :de_rand_1_bin_radiuslimited,
#         # DE/rand/2/bin:
#         :de_rand_2_bin,
#         # DE/rand/2/bin with radius limited sampling (a type of trivial geography):
#         :de_rand_2_bin_radiuslimited,
#         # Direct search:
#         # Generating set search:
#         # Compass/coordinate search:
#         :generating_set_search,
#         # Direct search through probabilistic descent:
#         :probabilistic_descent,
#         # Resampling Memetic Searchers:
#         # Resampling Memetic Search (RS):
#         :resampling_memetic_search,
#         # Resampling Inheritance Memetic Search (RIS):
#         :resampling_inheritance_memetic_search,
#         # Stochastic Approximation:
#         # Simultaneous Perturbation Stochastic Approximation (SPSA):
#         :simultaneous_perturbation_stochastic_approximation,
#     ]
#     results = []
#     for m in algs
#         println(m)
#         settings[:options][:Method] = m
#         _r = runfit(settings)
#         push!(results, _r)
#     end
#     results
# end

# begin
#     fitnesses = []
#     settings_file = joinpath(@__DIR__, "../data/default.jl")
#     include(settings_file)
#     settings[:options][:MaxTime] = 5
#     settings[:options][:TraceMode] = :silent
#     for _s in 1:4:41
#         println(_s)
#         settings[:options][:SamplerRadius] = _s
#         r = [runfit(settings) for _ in 1:5]
#         f = best_fitness.(r)
#         push!(fitnesses, sum(f))
#     end
# end
#
# begin
#     settings_file = joinpath(@__DIR__, "../data/default.jl")
#     gendata = joinpath(@__DIR__, "../test/gendata.jl")
#     include(settings_file)
#     include(gendata)
#
#     data = get_data(settings)
#     p = [A_ssp..., A_ppp..., ω..., Γ..., a_ssp'..., Δω..., a_pow..., χ3..., φ...]
#     p[1:4] ./= 10
#     p[5:6] ./= 10000
#     p[7:8] ./= 10
#     for i = 1:length(data)
#         y = model(data[i].ω, p, pol=data[i].pol, tstep=data[i].tstep, diff=data[i].diff)
#         @show i
#         @show data[i].signal[100], y[100]
#         @show data[i].signal ≈ y
#     end
# end
