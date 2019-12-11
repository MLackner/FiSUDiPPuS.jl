module FiSUDiPPuS

using Optim: optimize, LBFGS, ParticleSwarm, Fminbox
import Optim
using Printf: @printf
using Images: imresize
using JLD2
using FileIO: load, save
using Dates: now

export runfit, printresult, model, get_data

function runfit(options::String=joinpath(@__DIR__, "../data/default.jl");
        share_a::Bool=true, saveprefix="", multithreading=true)
    include(options)
    runfit(settings; share_a=share_a, saveprefix=saveprefix, multithreading=true)
end

function runfit(settings::Dict; saveprefix="")
    ### CONTENTS OF OPTIONS
    ## OPTIMIZATION
    # iterations
    # n_particles
    # lower
    # upper
    # start
    # time_limit
    ## DATA
    # filepath_ssp
    # filepath_ppp
    # roi_ssp
    # roi_ppp
    # size_ssp
    # size_ppp
    # reference_ssp
    # reference_ppp
    # bleach_weight
    # datatype
    # N
    ## SAVING
    # savedir

    ## Convert to correct types
    start = settings[:datatype].(settings[:start])
    lower = settings[:datatype].(settings[:lower])
    upper = settings[:datatype].(settings[:upper])
    bleach_weight = settings[:datatype](settings[:bleach_weight])

    # Process the data
    spectra = get_data(settings)

    c = Array{Float64,1}(undef, length(spectra)) # start cost
    function cost(p)
        Threads.@threads for i in 1:length(spectra)
            c[i] = sum(
                (model(spectra[i].ω, p; pol=spectra[i].pol,
                        diff=spectra[i].diff, tstep=spectra[i].tstep) .-
                spectra[i].signal).^2
            )
            if spectra[i].diff == true
                c[i] *= settings[:bleach_weight]
            end
        end
        sum(c)
    end

    if settings[:method] == ParticleSwarm()
        method = ParticleSwarm(
            lower=lower,
            upper=upper,
            n_particles=settings[:n_particles],
        )
    else
        method = Fminbox(settings[:method])
    end

    result = optimize(
        cost, lower, upper, start,
        method,
        Optim.Options(
            iterations=settings[:iterations],
            show_trace=false,
            time_limit=settings[:time_limit],
            )
        )

    save_result(settings[:savedir], spectra, result, settings;
                saveprefix=saveprefix)

    result
end

function viewsettings(file::String)
    include(file)   # contains the settings
    viewsettings(settings)
end

function viewsettings(s::Dict)
    function draw_roi(ax, roi)
        y0 = roi[1][1] - 1 # the plot starts at 0
        x0 = roi[2][1] - 1
        y1 = roi[1][end]
        x1 = roi[2][end]

        rectx = [
            x0 x1 x1 x0
            x1 x1 x0 x0
        ]
        recty = [
            y0 y0 y1 y1
            y0 y1 y1 y0
        ]

        ax.plot(rectx, recty, "r", linewidth=3)
    end

    x, y = get_data(s)

    # load the signal data to view the ROIs
    data_ssp = load(s[:filepath_ssp])
    data_ppp = load(s[:filepath_ppp])
    z_ssp = data_ssp["signal"]
    z_ppp = data_ppp["signal"]

    # Calculate start values
    ystart = startdata(settings)

    f1 = figure()
    plot(y, label="data")
    plot(ystart, label="initial")
    legend()

    f2 = figure()
    ax = subplot(111)
    ax.set_title("SSP Data ROI")
    ax.pcolormesh(z_ssp)
    draw_roi(ax, s[:roi_ssp])

    f3 = figure()
    ax = subplot(111)
    ax.set_title("PPP Data ROI")
    ax.pcolormesh(z_ppp)
    draw_roi(ax, s[:roi_ppp])

    f1, f2, f3
end

function get_data(s::Dict)
    get_data(
        s[:filepath_ssp],
        s[:filepath_ppp],
        s[:roi_ssp],
        s[:roi_ppp],
        s[:size_ssp],
        s[:size_ppp],
        s[:reference_ssp],
        s[:reference_ppp],
        s[:bleach_weight],
        s[:datatype],
    )
end

struct Spectrum
    ω
    signal
    pol::Symbol
    diff::Bool
    tstep::Int
end


"""
Returns an array of spectra
"""
function get_data(
        filepath_ssp,
        filepath_ppp,
        roi_ssp,
        roi_ppp,
        size_ssp,
        size_ppp,
        reference_ssp,
        reference_ppp,
        bleach_weight,
        datatype,
        )

    data_ssp = load(filepath_ssp)
    data_ppp = load(filepath_ppp)

    # time steps (this ends up being simply a unit range from
    # 1 to the size of the resized data)
    x1_ssp = 1:size_ssp[1]
    x1_ppp = 1:size_ppp[1]
    # wavenumbers
    x2_ssp = data_ssp["wavenumber"][roi_ssp[2]]
    x2_ppp = data_ppp["wavenumber"][roi_ppp[2]];
    # resize wavenumbers
    x2_ssp = imresize(x2_ssp, size_ssp[2])
    x2_ppp = imresize(x2_ppp, size_ppp[2])
    # signal
    y_ssp = data_ssp["signal"][roi_ssp...]
    y_ppp = data_ppp["signal"][roi_ppp...]

    # order the spectra from low to high delay in the roi
    p_ssp = sortperm(data_ssp["dltime"][roi_ssp[1]])
    p_ppp = sortperm(data_ppp["dltime"][roi_ppp[1]])
    # permute
    y_ssp = y_ssp[p_ssp,:]
    y_ppp = y_ppp[p_ppp,:]

    # resize signal
    y_ssp = imresize(y_ssp, size_ssp)
    y_ppp = imresize(y_ppp, size_ppp)

    s = Spectrum[] # put all spectra in here
    for (y,pol) in zip([y_ssp, y_ppp], [:ssp, :ppp])
        for i = 1:size(y,1)
            if pol == :ssp
                ω = x2_ssp
                tstep = x1_ssp[i]
            elseif pol == :ppp
                ω = x2_ppp
                tstep = x1_ppp[i]
            end
            for isdiff = [true, false]
                if !isdiff
                    # plain spectrum
                    signal = y[i,:]
                else
                    # difference spectrum
                    if pol == :ssp
                        r = reference_ssp
                    else
                        r = reference_ppp
                    end
                    signal = y[i,:] .- y[r,:]
                end
                _s = Spectrum(
                    ω,
                    signal,
                    pol,
                    isdiff,
                    tstep,
                )
                push!(s, _s)
            end
        end
    end
    return s
end


function get_results(filepath)
    data = load(filepath)
    x = data["x"]
    y = data["y"]
    r = data["result"]
    settings = data["settings"]

    x, y, r, settings
end

function plotresult(filepath::String)
    x,y,r,s = get_results(filepath)
    T = s[:datatype]

    ystart = model(
        x, T.(s[:start]),
        s[:N],
        T(s[:bleach_weight])
    )

    ymin = model(
        x, r.minimizer,
        s[:N],
        T(s[:bleach_weight])
    )

    f = figure()
    plot(y, label="data")
    plot(ystart, label="start")
    plot(ymin, label="minimized")
    legend()
    f
end

function printresult(r::String)
    data = load(r)
    result = data["result"]
    printresult(result)
end

function printresult(r)
    # Calculate the number of resonances by checking how many big-number
    # parameters there are.
    N = length(findall(x -> x > 1000, r.minimizer))

    # Get all the parameters from the result
    Assp = r.minimizer[   1: N]
    Appp = r.minimizer[ N+1:2N]
    ω    = r.minimizer[2N+1:3N]
    Γ    = r.minimizer[3N+1:4N]
    Δω   = r.minimizer[end-3]
    a_pow= r.minimizer[end-2]
    χ3   = r.minimizer[end-1]
    φ    = r.minimizer[end]

    # Get the initial parameter
    Assp0 = r.initial_x[1:N]
    Appp0 = r.initial_x[N+1:2N]
    ω0    = r.initial_x[2N+1:3N]
    Γ0    = r.initial_x[3N+1:4N]
    Δω0   = r.initial_x[end-3]
    a_pow0=r.initial_x[end-2]
    χ30   = r.initial_x[end-1]
    φ0    = r.initial_x[end]

    n = ["", "Assp", "Appp", "ω", "Γ", "A₀ssp", "A₀ppp", "ω₀", "Γ₀"]
    @printf "Δω_ppp: %.3f (start: %.3f)\n" Δω Δω0
    @printf "a_pow:  %.3f (start: %.3f)\n" a_pow a_pow0
    @printf "χ3:     %.3f (start: %.3f)\n" χ3 χ30
    @printf "φ:      %.3f (start: %.3f)\n" φ φ0

    # this prints the header
    for i = 1:length(n)
        @printf "%8.8s" n[i]
    end
    print("\n")

    # this prints the content of the table
    for i = 1:N
        @printf "%8.4u %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" i Assp[i] Appp[i] ω[i] Γ[i] Assp0[i] Appp0[i] ω0[i] Γ0[i]
    end

    nothing
end

"""
sum frequency spectrum
Keyword arguments:
φ:  χ(3) phase angle
χ3: χ(3) strength
"""
function sfspec(x, A, ω, Γ; φ=0π, χ3=0.0)
    ycmplx = 0.0 + 0.0im
    for i = 1:size(A,1)
        ycmplx += A[i] / (x - ω[i] - 1im * Γ[i])
    end
    ycmplx += χ3 * exp(1im * φ)
    abs2(ycmplx)
end

"""
p vector:
* A_ssp
* A_ppp
* ω
* Γ
* Δω_ppp
* a_pow
* χ3
* φ
"""
function model(x, p; pol=:none, diff=false, tstep=1)
    pol ≠ :ssp && pol ≠ :ppp && error("pol $pol not defined")

    # Get number of resonances
    N = (length(p) - 4) ÷ 6

    # Decompse parameter array
    A = Array{Float64,1}(undef,N)
    ω = Array{Float64,1}(undef,N)
    if pol == :ssp
        A .= p[   1: N]
        ω .= p[2N+1:3N]
    elseif pol == :ppp
        A .= p[ N+1:2N]
        ω .= p[2N+1:3N] .+ p[end-3] # this is the shift
    end
    Γ     = p[3N+1:4N]

    Δω_ppp= p[end-3]
    a_pow = p[end-2]
    χ3    = p[end-1]
    φ     = p[end]

    α = p[(3+tstep)*N+1:(4+tstep)*N]
    pol == :ppp && (α .^= a_pow)
    α .*= A

    # Preallocate signal
    y = Array{Float64,1}(undef, length(x))

    for i in eachindex(y)
        y[i] = sfspec(x[i], α, ω, Γ; φ=φ, χ3=χ3)
    end
    if diff == true
        for i in eachindex(y)
            y[i] -= sfspec(x[i], A, ω, Γ; φ=φ, χ3=χ3)
        end
    end

    y
end

function save_result(savedir,spectra,r,settings; saveprefix="")
    method = typeof(r.method)
    !isdir(savedir) && mkdir(savedir)
    datestr = string(now())
    datestr = replace(datestr, ":" => ".")
    filename = saveprefix * "_" * datestr * ".jld2"
    savepath = joinpath(savedir, filename)

    d = Dict(
        "result" => r,
        "spectra" => spectra,
        "settings" => settings,
    )

    save(savepath, d)
    nothing
end

end # module
