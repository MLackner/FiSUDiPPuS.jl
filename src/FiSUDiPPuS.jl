module FiSUDiPPuS

using BlackBoxOptim
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
    bleach_weight = settings[:datatype](settings[:bleach_weight])

    # Process the data
    spectra = get_data(settings)

    # Calculate universal cost. This is returned as an array. We have to then
    # check if the :borg_moea algorithm was selected or any other. The borg_moea
    # requires a tuple of fitnesses all other algorithms require a single float
    # value.
    function cost_universal(p)
        c = Array{Float64,1}(undef, length(spectra)) # start cost
        Threads.@threads for i in 1:length(spectra)
            c[i] = sum(
                (model(spectra[i].ω, p;
                        pol=spectra[i].pol,
                        diff=spectra[i].diff,
                        tstep=spectra[i].tstep,
                        N=settings[:N],
                        n_steps=settings[:n_steps],
                        ω_shift=settings[:ω_shift]) .-
                spectra[i].signal).^2
            )
            if spectra[i].diff == true
                c[i] *= settings[:bleach_weight]
            end
        end
        c
    end

    if settings[:options][:Method] == :borg_moea
        cost = p -> tuple(cost_universal(p)...)
    else
        cost = p -> sum(cost_universal(p))
    end

    result = bboptimize(
        cost;
        settings[:options]..., # solver options (bboptimize)
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

    # In case we have only one spectrum there is no difference spectrum provided
    # So we can delete the first and third spectra
    # if size_ssp[1] == 1 && size_ppp[1] == 1
    #     deleteat!(s, [1,3])
    # end
    # @show size(s)
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

function printresult(r::String, N)
    data = load(r)
    result = data["result"]
    printresult(result, N)
end

function printresult(r, N, ω_shift=false)
    p = best_candidate(r)

    # Get all the parameters from the result
    Assp = p[   1: N]
    Appp = p[ N+1:2N]
    ω    = p[2N+1:3N]
    Γ    = p[3N+1:4N]
    φssp = p[4N+1:5N]
    φppp = p[5N+1:6N]
    Δω   = p[end-2]
    a_pow= p[end-1]
    χnr  = p[end]

    # all step dependent parameters
    sdp = p[6N+1:end-3]
    n_sdp = length(sdp)
    ω_shift ? (n_steps = n_sdp ÷ 2 ÷ N) : (n_steps = n_sdp ÷ N)
    if ω_shift == true
        a  = sdp[1:n_steps*N]
        δω = sdp[n_steps*N+1:end]
    else
        a  = sdp
        δω = fill(0.0, n_sdp)
    end

    a  = reshape(a , (N, n_steps))'
    δω = reshape(δω, (N, n_steps))'

    n = ["", "Assp", "Appp", "ω", "Γ", "φssp", "φppp"]
    @printf "Δω_ppp: %.3f\n" Δω
    @printf "a_pow:  %.3f\n" a_pow
    @printf "χnr:    %.3f\n" χnr

    # this prints the header
    for i = 1:length(n)
        @printf "%9s" n[i]
    end
    print("\n")

    # this prints the content of the table
    for i = 1:N
        @printf "%8.4u %8.4f %8.4f %8.4f %8.4f %8.4fπ %8.4fπ\n" i Assp[i] Appp[i] ω[i] Γ[i] φssp[i]/π φppp[i]/π
    end

    print("\n")
    # print header for time dependet parameters
    @printf "%8.8s" "STEP"
    for i = 1:N
        header_a  = "a$i"
        header_δω = "δω$i"
        @printf "%8.8s %8.8s" header_a header_δω
    end
    print("\n")

    for i in 1:n_steps
        @printf "%8.4u" i
        for j in 1:N
            @printf "%8.4f %8.4f" a[i,j] δω[i,j]
        end
        print("\n")
    end

    nothing
end

"""
sum frequency spectrum
Keyword arguments:
χnr: nonresonant background
"""
function sfspec(x, A, ω, Γ, φ=zeros(length(A)); χnr=0.0)
    ycmplx = 0.0 + 0.0im
    for i = 1:size(A,1)
        ycmplx += A[i] / (x - ω[i] - 1im * Γ[i]) * exp(1im * φ[i])
    end
    ycmplx += χnr
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
function model(x, p; pol=:none, diff=false, tstep=1, N=1, ω_shift=false, n_steps=2)
    pol ≠ :ssp && pol ≠ :ppp && error("pol $pol not defined")

    # 01: A_ssp_1
    # 02: A_ssp_2
    # 03: A_ppp_1
    # 04: A_ppp_2
    # 05: ω_1
    # 06: ω_2
    # 07: Γ_1
    # 08: Γ_2
    # 09: φ_1
    # 10: φ_2
    # 11: a_11
    # 12: a_21
    # 13: a_12
    # 14: a_22
    # 15: Δω_11
    # 16: Δω_21
    # 17: Δω_12
    # 18: Δω_12
    # 19: Δω_ppp
    # 20: a_pow
    # 21: χnr

    # Decompse parameter array
    A = Array{Float64,1}(undef,N)
    ω = Array{Float64,1}(undef,N)
    φ = Array{Float64,1}(undef,N)
    if pol == :ssp
        A .= p[   1: N] .* 10
        ω .= p[2N+1:3N] .* 1e4
        φ .= p[4N+1:5N]
    elseif pol == :ppp
        A .= p[ N+1:2N] .* 10
        ω .= p[2N+1:3N] .* 1e4 .+ p[end-2] # this is the shift
        φ .= p[5N+1:6N]
    end
    Γ     = p[3N+1:4N] .* 10

    if ω_shift == true
        # We have extra parameters for the shift of each resonance
        δω = p[6N+N*n_steps+1+N*(tstep-1) : 7N+N*n_steps+N*(tstep-1)]
        ω .+= δω # apply wavelength shift
    end

    Δω_ppp= p[end-2]
    a_pow = p[end-1]
    χnr    = p[end]

    α = p[6N+1+N*(tstep - 1) : 7N+N*(tstep - 1)]
    pol == :ppp && (α .^= a_pow) # adjust for different pump strength
    α .*= A

    # Preallocate signal
    y = Array{Float64,1}(undef, length(x))

    for i in eachindex(y)
        y[i] = sfspec(x[i], α, ω, Γ, φ; χnr=χnr)
    end

    ω_shift && (ω .-= δω) # disapply wavelength shift because we use the
                          # reference spectrum in the next step
    if diff == true
        for i in eachindex(y)
            y[i] -= sfspec(x[i], A, ω, Γ, φ; χnr=χnr)
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
    :nothing
end

end # module
