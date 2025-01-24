""" run_trajectories(sys::System, params::SimulParameters) -> Vector{Trajectory}

Sample multiple trajectories for a given system and parameters.

# Arguments
- `sys::System`: The quantum system to simulate, containing information about its structure, energy levels, and dynamics.
- `params::SimulParameters`: A structure containing simulation parameters such as:
- `progbar::Bool`: show progress bar or not. `true` by default.

# Returns
- `Vector{Trajectory}`: A vector containing the results of the simulated trajectories. Each element corresponds to a single trajectory and encapsulates relevant system state information over time.
"""
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
            data[k] = run_single_trajectory(sys, params, W[:, tid], P[:, tid],
                        ts, Qs, Vs; seed = params.seed + k, isrenewal=isrenewal)
            next!(p)
            end
    finish!(p)
        return data
    # end
end

