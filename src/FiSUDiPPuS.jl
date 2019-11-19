module FiSUDiPPuS

using Optim: optimize, LBFGS, ParticleSwarm, Fminbox
import Optim
using Printf: @printf
using Images: imresize
using JLD2
using FileIO: load, save
using Dates: now
using PyPlot: figure, plot, legend, subplot

export runfit, viewsettings, plotresult, printresult

function runfit(options::String=joinpath(@__DIR__, "../data/default.jl"); share_a::Bool=true)
    include(options)
    runfit(settings)
end

function runfit(settings::Dict; share_a::Bool=true)
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
    x, y = get_data(
        settings[:filepath_ssp],
        settings[:filepath_ppp],
        settings[:roi_ssp],
        settings[:roi_ppp],
        settings[:size_ssp],
        settings[:size_ppp],
        settings[:reference_ssp],
        settings[:reference_ppp,],
        bleach_weight,
        settings[:datatype],
    )

    costfun(p) = sum((model(x, p, settings[:N], bleach_weight; share_a=share_a) .- y).^2)

    if settings[:method] == :ParticleSwarm
        method = ParticleSwarm(
            lower=lower,
            upper=upper,
            n_particles=settings[:n_particles],
        )
    elseif settings[:method] == :LBFGS
        method = Fminbox(LBFGS())
    end

    result = optimize(
        costfun, lower, upper, start,
        method,
        Optim.Options(
            iterations=settings[:iterations],
            show_trace=false,
            time_limit=settings[:time_limit],
            )
        )

    save_result(settings[:savedir], x, y, result, settings)

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

function startdata(settings::Dict)
    x, _ = get_data(settings)
    ystart = FiSUDiPPuS.model(
        x, settings[:datatype].(settings[:start]),
        settings[:N],
        settings[:datatype](settings[:bleach_weight])
        )
    ystart
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

    x = build_independents(x1_ssp, x1_ppp, x2_ssp, x2_ppp)
    y = build_dependents(y_ssp, y_ppp, reference_ssp, reference_ppp, bleach_weight)

    # println("size ssp data: $(size(data_ssp["signal"])) | region of interest: $roi_ssp | resized to $size_ssp")
    # println("size ppp data: $(size(data_ppp["signal"])) | region of interest: $roi_ppp | resized to $size_ppp")

    datatype.(x), datatype.(y)
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
    a_pow= r.minimizer[end-1]
    Δω   = r.minimizer[end]

    # Get the initial parameter
    Assp0 = r.initial_x[1:N]
    Appp0 = r.initial_x[N+1:2N]
    ω0    = r.initial_x[2N+1:3N]
    Γ0    = r.initial_x[3N+1:4N]
    a_pow0= r.initial_x[end-1]
    Δω0   = r.initial_x[end]

    n = ["", "Assp", "Appp", "ω", "Γ", "A₀ssp", "A₀ppp", "ω₀", "Γ₀"]
    @printf "Δω_ppp: %.3f (start: %.3f)\n" Δω Δω0
    @printf "a_pow:  %.3f (start: %.3f)\n" a_pow a_pow0

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
"""
function sfspec(x, A, ω, Γ)
    y = 0.0 + 0.0im
    for i = 1:size(A,1)
        y += A[i] / (x - ω[i] - 1im * Γ[i]) #* exp(1im * φ[i])
    end
    abs2(y)
end


function model(x::Array{T,2}, p::Array{T,1}, N::Int, bleach_weight::T;
               share_a::Bool=true) where T

    t = @view x[:,1]
    wn = @view x[:,2]
    pol = @view x[:,3]
    dif = @view x[:,4]

    A_ssp = p[   1: N]
    A_ppp = p[ N+1:2N]
    ω     = p[2N+1:3N]
    Γ     = p[3N+1:4N]

    a_pow = p[end-1]    # power by which the a's for the ppp spectra differ from
                        # those of the ssp spectra
    ω_ppp = ω .+ p[end] # add offset parameter to resonance positions for
                        # ppp spectra

    a_idx = 4N+1

    y = Array{T,1}(undef, size(x,1))
    α = Array{T,1}(undef, N)
    a = p[a_idx:a_idx+N-1]
    # define our α
    # TODO: Make sure the ssp spectrum comes first!
    α = a .* A_ssp

    @inbounds for i ∈ eachindex(y)
        if pol[i] == 0
            #ssp spectrum
            # α .= a .* A_ssp

            if dif[i] == 1
                # difference ssp
                y[i] = sfspec(wn[i], α, ω, Γ)
                y[i] -= sfspec(wn[i], A_ssp, ω, Γ)
                y[i] *= bleach_weight
            else
                # plain spectrum ssp
                y[i] = sfspec(wn[i], α, ω, Γ)
            end
        else
            # ppp spectrum
            # if share_a
            #     @. α = (a ^ a_pow) * A_ppp
            # else
            #     @. α = a * A_ppp
            # end

            if dif[i] == 1
                # difference ppp
                y[i] = sfspec(wn[i], α, ω_ppp, Γ)
                y[i] -= sfspec(wn[i], A_ppp, ω_ppp, Γ)
                y[i] *= bleach_weight
            else
                # plain spectrum ppp
                y[i] = sfspec(wn[i], α, ω_ppp, Γ)
            end
        end
        if i != length(y)
            # check if we have to increase the index for
            # the `a`s
            if share_a && t[i+1] < t[i]
                # reset to start value if we change to
                # the next spectrum (no matter if its a plain or
                # a difference spectrum)
                # TODO: make sure the spectra are ordered
                # correctly timewise
                a_idx = 4N+1
                a .= p[a_idx:a_idx+N-1]
                # Since the α values are dependent on the a values we
                # have to change them too. Additionaly we have to check if
                # there's an ssp or a ppp spectrum coming up next.
                if pol[i+1] == 0
                    @. α = a * A_ssp
                else
                    @. α = a ^ a_pow * A_ppp
                end
            elseif t[i+1] != t[i] # just go to the next time step
                a_idx += N
                a .= p[a_idx:a_idx+N-1]
                if pol[i+1] == 0
                    @. α = a * A_ssp
                else
                    @. α = a ^ a_pow * A_ppp
                end
            end
            if !share_a && dif[i+1] > dif[i]
                # if the a values shall not be shared between the
                # spectra reset to the start index when we switch
                # to the difference spectra
                a_idx = 4N+1
                a .= p[a_idx:a_idx+N-1]
                if pol[i+1] == 0
                    @. α = a * A_ssp
                else
                    @. α = a ^ a_pow * A_ppp
                end
            end
        end
    end
    y
end

function convert_independents(x1, x2)
    a1 = eltype(x1)[]
    a2 = eltype(x2)[]
    for _x1 in x1, _x2 in x2
        push!(a1, _x1)
        push!(a2, _x2)
    end
    x = [a1 a2];
end

function build_independents(x1_ssp, x1_ppp, x2_ssp, x2_ppp)
    # we have the following order top to bottom:
    #    all ssp spectra
    #    all ppp spectra
    #    all ssp difference spectra
    #    all ppp difference spectra

    # the convert_independents function gives just a set of
    # time steps and wavenumbers for each polarization
    x_ssp = convert_independents(x1_ssp, x2_ssp)
    x_ppp = convert_independents(x1_ppp, x2_ppp)

    # we can add the polarization information and that we have spectra
    x_ssp = [x_ssp fill(0.0, size(x_ssp,1))]
    x_ppp = [x_ppp fill(1.0, size(x_ppp,1))]

    # to these we can now add if it's a plain spectrum or if we have
    # a difference spectrum
    x_ssp_spec = [x_ssp fill(0.0, size(x_ssp,1))]
    x_ppp_spec = [x_ppp fill(0.0, size(x_ppp,1))]
    x_ssp_diff = [x_ssp fill(1.0, size(x_ssp,1))]
    x_ppp_diff = [x_ppp fill(1.0, size(x_ppp,1))]

    # now we just have to stich everything together in the correct order
    # (vertical concatenation)
    [x_ssp_spec; x_ppp_spec; x_ssp_diff; x_ppp_diff]
end

function build_dependents(y_ssp, y_ppp, r_ssp::Int, r_ppp::Int, bleach_weight)
    # first the plain spectra (turn the matrix and reshape to 1d array)
    # then concatenate vertically
    y_spec = [y_ssp'[:]; y_ppp'[:]]
    # get the reference spectra
    y_ref_ssp = y_ssp[r_ssp,:]
    y_ref_ppp = y_ppp[r_ppp,:]
    # calculate the difference spectra and concatenate
    y_diff_ssp = vcat([y_ssp[i,:] .- y_ref_ssp for i=1:size(y_ssp,1)]...) .* bleach_weight
    y_diff_ppp = vcat([y_ppp[i,:] .- y_ref_ppp for i=1:size(y_ppp,1)]...) .* bleach_weight
    y_diff = [y_diff_ssp; y_diff_ppp]
    [y_spec; y_diff]
end

function save_result(savedir,x,y,r,settings)
    method = typeof(r.method)
    !isdir(savedir) && mkdir(savedir)
    filename = string(now()) * ".jld2"
    savepath = joinpath(savedir, filename)

    d = Dict(
        "result" => r,
        "x" => x,
        "y" => y,
        "settings" => settings,
    )

    save(savepath, d)
    nothing
end

end # module
