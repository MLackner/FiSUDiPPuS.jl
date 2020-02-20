## DEFAULT SETTING FOR TESTING

# Data Settings
settings = Dict(
    :filepath_ssp => joinpath(@__DIR__, "sample_ssp.jld2"),
    :filepath_ppp => joinpath(@__DIR__, "sample_ppp.jld2"),
    # :roi_ssp => [1:3,1:301],        # region of interst [dl, wn]
    # :roi_ppp => [1:3,1:301],
    :roi_ssp => [1:2,1:301],
    :roi_ppp => [1:2,1:301],
    :size_ssp => (2,301),           # size of final data (dl, wn)
    :size_ppp => (2,301),
    :reference_ssp => 1,            # index of delays which's spectrum
    :reference_ppp => 1,            # should serve as reference
    :bleach_weight => 5,
    :datatype => Float64,
    :N => 3,                        # number of resonances

    :options => Dict(
        :TraceMode => :verbose,
        :TraceInterval => 1.0,
        :MaxSteps => 3e5,
        # :MaxTime => 20,
        :Method => :adaptive_de_rand_1_bin_radiuslimited,
        # :Method => :borg_moea,
        # :FitnessScheme => ParetoFitnessScheme{8}(is_minimizing=true),
        :MaxStepsWithoutProgress => 3e5,
        :ϵ => 0.101,
        :τ => 0.001, # above 0.001 is crap
        :θ => 0.90,
        :γ => 100.0, # no impact?
        :ζ => 1.0, # no impact?
        :RestartCheckPeriod => 10, # gets worse if higher
        :OperatorsUpdatePeriod => 10, # no impact?
        :NumRepetitions => 5, # doesnt work?
        :FitnessTolerance => 1e-11, # higher not good
        #:SamplerRadius => 8,
        :SearchRange => Tuple{Float64,Float64}[
            (-1, 1), # A1 ssp
            (-1, 1), # A2 ssp
            (-1, 1), # A3 ssp
            (-1, 1), # A1 ppp
            (-1, 1), # A2 ppp
            (-1, 1), # A3 ppp
            (0.2850, 0.2900), # ω
            (0.2900, 0.2920), # ω
            (0.2920, 0.2940), # ω
            (0.7, 1.2), # Γ
            (0.7, 1.2), # Γ
            (0.7, 1.2), # Γ
            fill((1.0, 1.0), 3)..., # a
            fill((0.8, 1.05), 3)..., # a
            (-5, 5), # Δω_ppp
            (0.1, 3.0), # a_pow
            (0, 0.4), # χ(3)
            (-2π, 2π), # φ
        ],
    ),
    :savedir => joinpath(@__DIR__, "results"), # saveing
)
