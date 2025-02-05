"""

```
run_trajectories(sys::System, params::SimulParameters; progbar::Bool = true,
                          psireset::VecOrMat{ComplexF64}=zeros(ComplexF64, 0))
```

Sample multiple trajectories of `sys`. with the `seed`, `nsamples` and `psi0` specified in
`params`. If the process is of renewal type,i.e. after a jump it always comes to the same
state, that state can be specified via `psireset` and it will be used to optimize the
sampling of the jump times.
#  Arguments
- `sys::System`: The system for which to run the trajectoris
- `params::SimulParameters`: specificies the details of the trajectiories:
                             seed, number of points in the finegrid, initial state
                             and tolerance for the dark state detection, and the number
                             of trajectories to run.

# Optional Arguments
- `progbar::Bool`: if true, show a progress bar of the iteration, `true` by default.
- `psireset::VecOrMat{ComplexF64}`: if specified, it is used to optimize the jump time sampling.

# Return
A `Vector{Trajectory}` of length `params.ntraj` with the sampled trajectories.
"""
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

