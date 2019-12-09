## DEFAULT SETTING FOR TESTING

# Data Settings
settings = Dict(
    :filepath_ssp => joinpath(@__DIR__, "sample_ssp.jld2"),
    :filepath_ppp => joinpath(@__DIR__, "sample_ppp.jld2"),
    # :roi_ssp => [1:3,1:301],        # region of interst [dl, wn]
    # :roi_ppp => [1:3,1:301],
    :roi_ssp => [1:2,10:290],
    :roi_ppp => [1:2,10:290],
    :size_ssp => (2,281),           # size of final data (dl, wn)
    :size_ppp => (2,281),
    :reference_ssp => 1,            # index of delays which's spectrum
    :reference_ppp => 1,            # should serve as reference
    :bleach_weight => 3,
    :datatype => Float64,
    :N => 3,                        # number of resonances

    # Optimization Settings
    :method => :BFGS,
    :iterations => 100_000_000_000,
    :n_particles => 50,
    :start => [
        -3, 4, -9,
        7, -5, 4,
        2872, 2911, 2938,
        fill(8.0, 3)...,
        fill(0.99, 3*2)...,
        1.0,        # Δω_ppp
        0.95,       # a_pow
        0.2,        # χ3
        -π/2        # φ
    ],
    :lower => [
        fill(-10, 3)...,
        fill(-10, 3)...,
        fill(2860, 3)...,
        fill(6.0, 3)...,
        fill(0.9, 3*2)...,
        0.0,
        -5.0,
        0.0,
        -4π,
    ],
    :upper => [
        fill(10, 3)...,
        fill(10, 3)...,
        fill(2950, 3)...,
        fill(12.0, 3)...,
        fill(1.05, 3*2)...,
        3.0,
        5.0,
        1.0,
        2π,
    ],
    :time_limit => 360,

    # Saving
    :savedir => joinpath(@__DIR__, "results"),
)
