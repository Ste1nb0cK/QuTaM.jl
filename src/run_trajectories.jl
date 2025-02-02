function run_trajectories(sys::System, params::SimulParameters; progbar::Bool = true, isrenewal::Bool=false)
    ## Precomputing
    t0 = params.dt # To avoid having jumps at 0
    ts = collect(LinRange(t0, params.multiplier*params.tf, params.nsamples))
    Qs = Array{ComplexF64}(undef, sys.NLEVELS, sys.NLEVELS, params.nsamples)
    Vs = Array{ComplexF64}(undef, sys.NLEVELS, sys.NLEVELS, params.nsamples)
    precompute!(sys, params.nsamples, ts, Qs, Vs)
    # To store the data
    data = Vector{Trajectory}(undef, params.ntraj)
    # Create where to store the weights, one for each thread.
    # In case the process is renewal a single copy of W is used
    W = Array{Float64}(undef, params.nsamples, isrenewal ? 1 : nthreads())
    P = Array{Float64}(undef, sys.NCHANNELS, nthreads())

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

