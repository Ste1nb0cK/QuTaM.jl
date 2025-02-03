function run_trajectories(sys::System, params::SimulParameters; progbar::Bool = true,
                          psireset::VecOrMat{ComplexF64}=zeros(ComplexF64, 0))
    ## Precomputing
    t0 = params.dt # To avoid having jumps at 0
    ts = collect(LinRange(t0, params.multiplier*params.tf, params.nsamples))
    Qs = Array{ComplexF64}(undef, sys.NLEVELS, sys.NLEVELS, params.nsamples)
    Vs = Array{ComplexF64}(undef, sys.NLEVELS, sys.NLEVELS, params.nsamples)
    precompute!(sys, params.nsamples, ts, Qs, Vs)
    # To store the data
    data = Vector{Trajectory}(undef, params.ntraj)
    # If a psireset was passed, it is assumed the process is of the renewal type
    isrenewal = !isempty(psireset)
    # Create where to store the weights, one for each thread.
    P = Array{Float64}(undef, sys.NCHANNELS, nthreads())
    # in case the process is renewal, call run_singletrajectory_renewal to optimize
    if isrenewal
        # Calculate the WTD weights of the initial state
        W0 = Vector{Float64}(undef, params.nsamples)
        calculatewtdweights!(W0, Qs, params.psi0, params)
        # and those of the reset state
        W = Vector{Float64}(undef, params.nsamples)
        calculatewtdweights!(W, Qs, psireset, params)
        # To avoid race conditions, pass to each thread a copy of the psireset
        psireset_copies = repeat(psireset, 1, nthreads())
        # now run the trajectories
        p = Progress(params.ntraj; dt=1.0, desc="Sampling...", enabled=progbar, showspeed=true)
        if ndims(psireset) == 1
            @threads for k in 1:params.ntraj
                tid  = threadid()
                data[k] = run_singletrajectory_renewal(sys, params,
                                           W, W0,
                                           P[:, tid],
                        ts, Qs, Vs, psireset_copies[:, tid]; seed = params.seed + k)
                next!(p)
            end
            finish!(p)
            return data
        elseif ndims(psireset == 2)
            @threads for k in 1:params.ntraj
                tid  = threadid()
                data[k] = run_singletrajectory_renewal(sys, params,
                                           W, W0,
                                           P[:, tid],
                        ts, Qs, Vs, psireset_copies[:, :, tid]; seed = params.seed + k)
                next!(p)
            end
            finish!(p)
            return data

        end


    end
    W = Array{Float64}(undef, params.nsamples,  nthreads())
    p = Progress(params.ntraj; dt=1.0, desc="Sampling...", enabled=progbar, showspeed=true)
   @threads for k in 1:params.ntraj
            tid  = threadid()
            data[k] = run_singletrajectory(sys, params,
                                           W[:, isrenewal ? 1 : tid],
                                           P[:, tid],
                        ts, Qs, Vs; seed = params.seed + k)
            next!(p)
            end
    finish!(p)
        return data
    # end
end

