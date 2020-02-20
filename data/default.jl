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
    :n_steps => 2,
    :datatype => Float64,
    :N => 3,                        # number of resonances
    # should the problem be solved allowing the resonances to shift?
    :ω_shift => true,

    :options => Dict(
        :TraceMode => :compact,
        :TraceInterval => 4.0,
        :MaxSteps => 3e7,
        :MaxStepsWithoutProgress => 3e7,
        :FitnessTolerance => 1e-11, # higher not good
        :PopulationSize => 1000, # :MaxTime => 20,

        # BORG MOEA
        # :Method => :borg_moea,
        # :FitnessScheme => ParetoFitnessScheme{8}(is_minimizing=true),
        :ϵ => 0.1,
        :τ => 0.001, # above 0.001 is crap
        :θ => 0.90,
        :γ => 4.0, # no impact?
        :ζ => 1.0, # no impact?
        :RestartCheckPeriod => 1000, # 10 is good gets worse if higher
        :OperatorsUpdatePeriod => 2, # no impact?
        :NumRepetitions => 5, # doesnt work?

        # SINGLE OBJECTIVE METHOD
        :Method => :adaptive_de_rand_1_bin_radiuslimited,
        :SearchRange => Tuple{Float64,Float64}[
            # A_ssp
            # (0.7, 1.1),
            # (-0.5, 0.0),
            # (-1.0, 1.0),
            # A_ppp
            # (0, 1.0),
            # (-1.0, 0.0),
            # (-1.0, 1.0),
            fill((-1.5, 1.5), 6)...,

            # (0.2875, 0.2980), # ω
            # (0.2880, 0.2905), # ω
            # (0.2920, 0.2940), # ω
            fill((0.2870, 0.2980), 3)..., # ω

            # Γ VALUES
            # (0.5, 2.0), # Γ
            # (0.3, 2.0), # Γ
            # (0.5, 2.0), # Γ
            fill((0.5, 2.0), 3)...,

            # a VALUES
            fill((1.0, 1.0), 3)..., # a
            fill((0.7, 1.15), 3)..., # a

            # δω VALUES
            fill((0.0, 0.0), 3)...,
            (0.0, 0.0),
            (-3.0, 1.0),
            (0.0, 0.0),

            # ADDITIONAL PARAMETERS
            (-5, 5), # Δω_ppp
            (0.1, 1.8), # a_pow
            (0, 0.1), # χ(3)
            (-2π, 2π), # φ
        ],
    ),
    :savedir => joinpath(@__DIR__, "results"), # saveing
)
