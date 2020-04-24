module FiSUDiPPuS

using BlackBoxOptim
using Printf: @sprintf
using Images: imresize
using JLD2
using FileIO: load, save
using Dates: now
using Statistics

export runfit, printresult, model, get_data, Spectrum, getparam


struct Spectrum
    ω
    signal
end


function runfit(spectra::Array{Spectrum}, settings::Dict)

    n_steps = settings[:n_steps]
    tsteps = Int[]
    specidx = Int[]
    for i in eachindex(spectra)
        push!(tsteps, 0:n_steps...)
        push!(specidx, fill(i, n_steps+1)...)
    end

    function cost(p)
        c = Array{Float64,1}(undef, length(spectra) + length(spectra) * n_steps) # start cost
        Threads.@threads for i in eachindex(c)
            c[i] = sum(
                (model(spectra[specidx[i]].ω, p;
                        specidx=specidx[i],
                        tstep=tsteps[i],
                        N=settings[:N],
                        n_steps=settings[:n_steps],
                        ω_shift=settings[:ω_shift]) .-
                spectra[specidx[i]].signal[tsteps[i]+1,:]).^2
            )
            if tsteps[i] > 0
                c[i] *= settings[:bleach_weight]
            end

        end
        sum(c)
    end

    result = bboptimize(
        cost;
        settings[:options]..., # solver options (bboptimize)
        )

    save_result(spectra, result, settings)
    printresult(result, settings, length(spectra)) |> print

    result
end



function get_results(filepath)
    data = load(filepath)
    x = data["x"]
    y = data["y"]
    r = data["result"]
    settings = data["settings"]

    x, y, r, settings
end

function printresult(result, settings, nspectra)
    @show p = best_candidate(result)

    # First we have to get the universal parameters
    # A, φ, ω, Γ, a, δω, χnr, Δω, β = decompose_parameters(p;
    x = [
        decompose_parameters(p;
            specidx=i,      # index of spectrum
            phase=settings[:phase],   # phase is not multiple of π?
            tstep=j,        # current time step
            N=settings[:N],            # number of resonances
            ω_shift=settings[:ω_shift],
            n_steps=settings[:n_steps],      # number of time steps
            ) for i in 1:nspectra, j in 0:settings[:n_steps]
        ]

    χnr = [x[i,1][7] for i in 1:nspectra]
    Δω  = [x[i,1][8] for i in 1:nspectra]
    β   = [x[i,settings[:n_steps]+1][9] for i in 1:nspectra]
    sp = [χnr..., Δω..., β...]

    str = ""

    universal_header = ["", "ω", "Γ"]
    per_spectrum_header = ["A", "φ"]
    spectrum_params_header = ["", "χₙᵣ", "Δω", "β"]
    steps_header = ["", "a", "δω"]

    # just that is there just once per spectrum
    for s in spectrum_params_header
        str *= @sprintf "%9s" s
    end

    str *= "\n"

    for i in 1:nspectra
        str *= @sprintf "%9.2u" i
        str *= @sprintf "%9.2f" χnr[i]
        str *= @sprintf "%9.2f" Δω[i]
        str *= @sprintf "%9.2f" β[i]

        str *= "\n"
    end

    str *= "\n"

    # this prints the header
    for s in universal_header
        str *= @sprintf "%9s" s
    end

    for i in 1:nspectra
        str *= @sprintf "%9s" "A" * string(i)
    end
    if settings[:phase]
        for i in 1:nspectra
            str *= @sprintf "%9s" "φ" * string(i)
        end
    end

    str *= "\n"

    for i in 1:settings[:N]
        # print resonance number
        str *= @sprintf "%9.2u" i
        # ω ond Γ
        ω = x[1,1][3][i]
        Γ = x[1,1][4][i]
        str *= @sprintf "%9.1f" ω
        str *= @sprintf "%9.2f" Γ

        # print A
        for j in 1:nspectra
            A = x[j,1][1][i]
            str *= @sprintf "%9.2f" A
        end
        # print phases
        if settings[:phase]
            for j in 1:nspectra
                φ = x[j,1][2][i]
                str *= @sprintf "%8.2fπ" φ / π
            end
        end

        str *= "\n"
    end

    str *= "\n"

    str *= @sprintf "%9s" steps_header[1]
    for n in 1:settings[:N]
        str *= @sprintf "%9s" steps_header[2] * string(n)
    end

    if settings[:ω_shift]
        for n in 1:settings[:N]
            str *= @sprintf "%9s" steps_header[3] * string(n)
        end
    end

    str *= "\n"

    for i in 1:settings[:n_steps]+1
        # step number
        str *= @sprintf "%9.3u" i-1
        # a value
        for n in 1:settings[:N]
            a = x[1,i][5][n]
            str *= @sprintf "%9.3f" a
        end

        if settings[:ω_shift]
            for n in 1:settings[:N]
               δω = x[1,i][6][n]
                str *= @sprintf "%9.3f" δω
            end
        end
        str *= "\n"
    end

    return str
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

function decompose_parameters(p;
        specidx=1,      # index of spectrum
        phase=false,   # phase is not multiple of π?
        tstep=0,        # current time step
        N=1,            # number of resonances
        ω_shift=false,
        n_steps=0,      # number of time steps
    )

    ### First Spectrum ###
    # 01: A_ik
    # 02: A_il
    # 03: φ_ik
    # 04: φ_il
    # 13: χnr
    # 05: Δω_i
    # 06: β_i
    ### Second Spectrum ###
    # (see above)
    ### Universal Parameters ###
    # 07: ω_k
    # 08: Γ_k
    # 05: a_kp
    # 06: a_kq
    # 07: a_lp
    # 08: a_lq
    # 09: δω_kp
    # 10: δω_kq
    # 11: δω_lp
    # 12: δω_lq

    # relevant parameters for this spectrum
    # parameters per spectrum (not including universal ones)
    pps = N + phase*N + 3
    qs = (specidx-1) * pps + 1 # parameters start
    qe = specidx * pps
    q = p[qs:qe]
    # universal parameters
    u = p[end-2N - n_steps*N - ω_shift*n_steps*N+1:end]

    # Decompse parameter array
    A = Array{Float64,1}(undef,N)
    φ = Array{Float64,1}(undef,N)
    ω = Array{Float64,1}(undef,N)
    Γ = Array{Float64,1}(undef,N)
    a = Array{Float64,1}(undef,N)
    δω = Array{Float64,1}(undef,N)

    A .= q[1:N]
    ω .= u[1:N]
    Γ .= u[N+1:2N]
    χnr = q[end-2]
    Δω =  q[end-1]
    β  =  q[end]

    if phase == true
        φ .= q[N+1:2N]
    else
        φ .= 0.0
    end

    if tstep > 0
        s = 2N + N * (tstep - 1) + 1
        e = s - 1 + N
        a .= u[s:e]
        if ω_shift
            s = 2N + n_steps*N + N * (tstep - 1) + 1
            e = s - 1 + N
            δω .= u[s:e]
        else
            δω .= 0.0
        end
    else
        a .= 1.0
        δω .= 0.0
        β = 1.0
    end

    return A, φ, ω, Γ, a, δω, χnr, Δω, β
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
function model(x, p; tstep=0, kwargs...)

    A, φ, ω, Γ, a, δω, χnr, Δω, β = decompose_parameters(p; tstep=tstep, kwargs...)

    # Preallocate signal
    y = Array{Float64,1}(undef, length(x))
    for i in eachindex(y)
        y[i] = sfspec(x[i], (a.^β) .* A, ω .+ Δω .+ δω, Γ, φ; χnr=χnr)
    end

    if tstep > 0
        for i in eachindex(y)
            y[i] -= sfspec(x[i], A, ω .+ Δω, Γ, φ; χnr=χnr)
        end
    end

    y
end

function save_result(spectra, result, settings)
    datestr = string(now())[1:19]
    datestr = replace(datestr, ":" => ".")
    savedir = "data/results"
    filename = datestr * ".jld2"
    savepath = joinpath(savedir, filename)

    !isdir(savedir) && mkdir(savedir)

    d = Dict(
        "result" => result,
        "spectra" => spectra,
        "settings" => settings,
    )

    save(savepath, d)
    :nothing
end

"""
Load a specific parameter from the result
tstep either :all or Int
"""
function getparam(result, settings, param::Symbol; tstep=0, specidx=1)
    if tstep==:all
        tstep = 1:settings[:n_steps]
    end
    pk = [:A, :φ, :ω, :Γ, :a, :δω, :χnr, :Δω, :β]
    p = best_candidate(result)
    param_out = zeros(settings[:N], length(tstep))
    for (i,t) in enumerate(tstep)
        pd = FiSUDiPPuS.decompose_parameters(p;
            specidx=specidx,
            phase=settings[:phase],
            ω_shift=settings[:ω_shift],
            tstep=t,
            N=settings[:N],
            n_steps=settings[:n_steps],
        )
        param_select = pd[pk .== param][1]
        param_out[:,i] .= param_select
    end
    param_out
end
getparam(data, param; kwargs...) = getparam(
    data["result"], data["settings"], param; kwargs...
)


end # module
