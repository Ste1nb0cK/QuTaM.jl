############# Precomputation routine ######################
"""

    precompute!(sys::System, nsamples::Int64,
         ts::Vector{Float64}, Qs::Vector{Matrix{ComplexF64}})

Does the precomputation routine for the Gillipsie algorithm.
The result is stored in the `Qs`.
# Arguments

`sys::System`: system
`nsamples::Int64`: number of samples
`ts::Vector{Float64}`: fine grid vector. *IT IS ASSUMED HOMOGENOUS*
`Qs::Vector{Matrix{ComplexF64}}`: vector of matrices to store the precomputation.
"""
function precompute!(sys::System, nsamples::Int64,
         ts::Vector{Float64}, Qs::Vector{Matrix{ComplexF64}})
        for k in 1:nsamples
            expm = exp(-1im*ts[k]*sys.Heff)
            Qs[k] = expm * sys.J * adjoint(expm)
        end
    return
end
############# Single Trajectory Routine ######################
"""
    run_single_trajectories(sys::System, params::SimulParameters,
     W::Vector{Float64}, P::Vector{Float64}, psi::Vector{ComplexF64},
     ts::Vector{Float64}, Qs::Vector{Matrix{ComplexF64}};
      seed::Int64 = 1) -> Trajectory

Sample a single trajectory from the system and parameters using the Gillipsie algorithm
. This is inteded to be used by `run_trajectories`

# Arguments:
- `sys::System`: System of interest
- `params::SimulParameters`: simulation parameters
- `W::Vector{Float64}`: to store the weights over the fine grid. Its lenght must be 1 more than that of ts
- `P::Vector{Float64}`: to store the weights over the channels
- `psi::Vector{ComplexF64}`: to store the current state vector
- `ts::Vector{Float64}`: the finegrid of waiting times
- `Qs::Vector{Matrix{ComplexF64}}`: to store the precomputed values
- `seed::Int64`: seed for generating the trajectory
# Returns:
The sample trajectory

# Warning: the seed does not coincide with that of `params` by default.

# Waning: final jump in the trajectory happens after final time
The trajectory ends when a jump that happens after `params.tf` is obtained,
yet that jump is stored in the trajectory. In other words, the last jump of the
trajectory always happen after the set final time.
"""
function run_single_trajectory(
    sys::System,
    params::SimulParameters,
    W::Vector{Float64}, P::Vector{Float64}, psi::Vector{ComplexF64}, ts::Vector{Float64},
    Qs::Vector{Matrix{ComplexF64}}; seed::Int64 = 1)
    # Random number generator
    Random.seed!(seed)
    traj = Vector{DetectionClick}()
    psi .= params.psi0
    t::Float64 = 0
    channel = 0
    # Run the trajectory
    while t < params.tf
        # Calculate the probability at infinity
        for k in 1:params.nsamples
           W[k] = real(dot(psi, Qs[k]*psi))
        end
        if sum(W) < params.eps
            break
        end
        # 2. Sample jump time
        tau = StatsBase.sample(ts, StatsBase.weights(W))
        t = tau + t
        psi .= exp(-1im*tau*sys.Heff) * psi
        # 3. Sample the channel
        aux_P = real(dot(psi, sys.J * psi))
        for k in 1:sys.NCHANNELS
            P[k] = norm(sys.Ls[k]*psi)^2
        end
        P .= P / aux_P
        channel::Int64 = StatsBase.sample(1:sys.NCHANNELS, StatsBase.weights(P))
        psi .= sys.Ls[channel]*psi # State without normalization
        psi .= psi / norm(psi)
        push!(traj, DetectionClick(tau, channel))
    end
    return traj
end

"""

    sample_single_trajectory(sys::System, params::SimulParameters, seed::Int64) -> Vector{Trajectory}

Sample a single trajectory. This is inteded to be used in case a single sample is needed without
redefining `params.ntraj`.

# Arguments
- `sys::System`
- `params::SimulParameters`
- `seed::Int64`

# Returnrs
A 1-element trajectory vector.
"""
function sample_single_trajectory(sys::System, params::SimulParameters, seed::Int64)
    ## Precomputing
    ts = collect(LinRange(0, params.multiplier*params.tf, params.nsamples))
    Qs = Vector{Matrix{ComplexF64}}(undef, params.nsamples)
    precompute!(sys, params.nsamples, ts, Qs)
    # Running the trajectory
    psi = Vector{ComplexF64}(undef, sys.NLEVELS)
    W = Vector{Float64}(undef, params.nsamples)
    P = Vector{Float64}(undef, sys.NCHANNELS)
    data = run_single_trajectory(sys, params,
                                W, P, psi, ts, Qs, seed = seed)
    return data
end

"""
    run_trajectories(sys::System, params::SimulParameters) -> Vector{Trajectory}

Sample multiple trajectories for a given system and parameters.

# Arguments
- `sys::System`: The quantum system to simulate, containing information about its structure, energy levels, and dynamics.
- `params::SimulParameters`: A structure containing simulation parameters such as:
- `progbar::Bool`: show progress bar or not. `true` by default.

# Returns
- `Vector{Trajectory}`: A vector containing the results of the simulated trajectories. Each element corresponds to a single trajectory and encapsulates relevant system state information over time.
"""
function run_trajectories(sys::System, params::SimulParameters; progbar::Bool = true)
    ## Precomputing
    t0 = eps(Float64) # To avoid having jumps at 0
    ts = collect(LinRange(t0, params.multiplier*params.tf, params.nsamples))
    Qs = Vector{Matrix{ComplexF64}}(undef, params.nsamples)
    precompute!(sys, params.nsamples, ts, Qs)
    # To store the data
    data = Vector{Trajectory}(undef, params.ntraj)
    # Running the trajectory
    psi = Vector{ComplexF64}(undef, sys.NLEVELS)
    W = Vector{Float64}(undef, params.nsamples)
    P = Vector{Float64}(undef, sys.NCHANNELS)
    if progbar
        @showprogress 1 "Sampling..." for k in 1:params.ntraj
            data[k] = run_single_trajectory(sys, params,
                                            W, P, psi, ts, Qs, seed = params.seed + k)
        end
        return data

    else
    for k in 1:params.ntraj
            data[k] = run_single_trajectory(sys, params,
                                            W,
                                            P,
                                            psi,
                                            ts, Qs, seed = params.seed + k)
        end
        return data

    end
end

############ Evaluation at given times #######################
"""
    states_at_jumps(traj::Trajectory, sys::System,
                      psi0::Vector{ComplexF64}) - > Array{ComplexF64}
From a given trajectory, recover the states at each jump.

# Arguments
- `traj::Trajectory`: Trajectory
- `sys::System`: System of interest
- `psi0::Vector{ComplexF64}`: Initial state
# Optional Arguments
- `normalize::Bool`: Whether to normalize the states or not, true by default

# Returns
`Array{ComplexF64}` with the states

The dimensions of the returned array `s` are `(size(traj), sys.NLEVELS)`,
so to recover the state vector at the ``n``-th jump one would do `s[n, :]`.
"""
function states_at_jumps(traj::Trajectory, sys::System,
                      psi0::Vector{ComplexF64}; normalize::Bool=true)
    njumps = size(traj)[1]
    # states = Vector{Vector{ComplexF64}}(undef, njumps)
    states = Array{ComplexF64}(undef, njumps, sys.NLEVELS)
    psi = copy(psi0)
    jump_counter = 1
    if normalize
        for click in traj
            psi .= sys.Ls[click.label] * exp(-1im*(click.time)*sys.Heff) * psi
            psi .= psi/norm(psi)
            for n in 1:sys.NLEVELS
                states[jump_counter, n] = psi[n]
            end
            jump_counter = jump_counter + 1
        end
        return states
    else
        for click in traj
            psi .= sys.Ls[click.label] * exp(-1im*(click.time)*sys.Heff) * psi
            for n in 1:sys.NLEVELS
                states[jump_counter, n] = psi[n]
            end
            jump_counter = jump_counter + 1
        end
        return states

    end
end

"""

    evaluate_at_t(t_given::Vector{Float64}, traj::Trajectory, sys::System,
                       psi0::Vector{ComplexF64}) -> Array{ComplexF64}

Evaluate in between jumps of the given trajectory and initial state.
The returned states are stored in a `Array{ComplexF64}` with dimensions
(size(t_given), sys.NLEVELS).

# Arguments
- `t_given::Vector{Float64}`: times at which the trajectory is to be evalauted
- `traj::Trajectory`: the trajectory
- `sys::System`: the system to which the trajectory corresponds
- `psi0::Vector{ComplexF64}`: the initial state of the trajectory
- `normalize::Bool`: whether to normalize the states or not, true by default
# Returns
A complex two-dimensional array whose rows contain the states.
"""

function evaluate_at_t(t_given::Vector{Float64}, traj::Trajectory, sys::System,
                       psi0::Vector{ComplexF64}; normalize::Bool=true)
    psi = copy(psi0)
    ntimes = size(t_given)[1]
    jump_states = states_at_jumps(traj, sys, psi0; normalize=normalize)
    njumps = size(jump_states)[1]
    t_ = 0
    counter = 1
    counter_c = 1
    # Special case: if the time array is empty, return an empty array
    if isempty(t_given)
        return Array{ComplexF64}(undef, 0, 0) # empty 2 dimensional array
    end

    states = Array{ComplexF64}(undef, ntimes, sys.NLEVELS)

    # Edge case: if the trajectory is empty, evaluate exponentials and return
    if isempty(traj)
        while counter <= ntimes
            psi .= exp(-1im*(t_given[counter])*sys.Heff) * psi
            if normalize
                psi .= psi/norm(psi)
            end
            for k in 1:sys.NLEVELS
                states[counter, k] = psi[k]
            end
            counter = counter + 1
            if counter > ntimes
                break
            end
        end
        return states
    end
    # All the states before the first jump can be handled like this:
    while (t_given[counter] < traj[counter_c].time) && (counter <= ntimes)
            psi .= exp(-1im*(t_given[counter])*sys.Heff) * psi0
            if normalize
                psi .= psi/norm(psi)
            end
            for k in 1:sys.NLEVELS
                states[counter, k] = psi[k]
            end
            counter = counter + 1
            if counter > ntimes
                break
            end
    end
    t_ = t_ + traj[counter_c].time
    counter_c = counter_c + 1
    while (counter_c <= njumps) && (counter <= ntimes)
        timeclick = traj[counter_c].time
        while (t_ < t_given[counter] < t_ + timeclick) && (counter <= ntimes)
             psi .= exp(-1im*(t_given[counter] - t_)*sys.Heff) * jump_states[counter_c-1, :]
             if normalize
             psi .= psi/norm(psi)
             end
             for k in 1:sys.NLEVELS
                 states[counter, k] = psi[k]
             end
             counter = counter + 1
             if counter > ntimes
                 break
             end
         end
       t_ = t_ + timeclick
       counter_c = counter_c + 1
    end

    while counter <= ntimes
        psi .= exp(-1im*(t_given[counter] - t_)*sys.Heff) * jump_states[njumps, :]
        if normalize
            psi .= psi/norm(psi)
        end
        for k in 1:sys.NLEVELS
            states[counter, k] = psi[k]
        end
        counter = counter + 1
    end
    return states
end


function DerivativeAtJumps(traj::Trajectory, sys::System, dHe::Matrix{ComplexF64}, dLs::Vector{Matrix{ComplexF64}},
         params::SimulParameters)
    njumps = size(traj)[1]
    dpsis = zeros(ComplexF64, njumps, sys.NLEVELS)
    tmp1 = zeros(ComplexF64, sys.NLEVELS) #for the psitilde
    tmp2 = zeros(ComplexF64, sys.NLEVELS) #for the derivative of psitilde
    #  intialize the monitoring
    label = traj[1].label
    tau = traj[1].time
    tmp1 .= params.psi0
    #  Iterate over the clicks
    for k in 1:njumps
        label = traj[k].label
        tau = traj[k].time
        # 1. Calculate the derivative
        tmp2 .= dLs[label]*exp(-1im*tau*sys.Heff)*tmp1 +
                sys.Ls[label]*-1im*tau*dHe*exp(-1im*tau*sys.Heff)*tmp1+
                sys.Ls[label]*exp(-1im*tau*sys.Heff)*tmp2
        dpsis[k, :] = tmp2
        # 2. Calculate the psitilde
        tmp1 .= sys.Ls[label]*exp(-1im*tau*sys.Heff)*tmp1
    end
    return dpsis
end


function MonitoringInBetween(
        traj::Trajectory, sys::System, dHe::Matrix{ComplexF64}, dLs::Vector{Matrix{ComplexF64}},
         params::SimulParameters, t_given::Vector{Float64})
    ntimes = size(t_given)[1]
    tmp1 = zeros(ComplexF64, sys.NLEVELS)
    psitilde = zeros(ComplexF64, sys.NLEVELS)
    # Obtain the unnormalized jump states
    jump_states = states_at_jumps(traj, sys, params.psi0; normalize=false)
    inbetween_states = evaluate_at_t(t_given, traj, sys, params.psi0; normalize=false)
    jump_dpsis= DerivativeAtJumps(traj, sys, dHe, dLs, params)
    njumps = size(jump_states)[1]
    t_ = 0
    counter = 1
    counter_c = 1
    # Edge case: if the time array is empty, return an empty array
    if isempty(t_given)
        return Array{ComplexF64}(undef, 0, 0) # empty 2 dimensional array
    end

    xis = Array{ComplexF64}(undef, ntimes, sys.NLEVELS, sys.NLEVELS)
    #Edge case: if the trajectory is empty
     if isempty(traj)
        while counter <= ntimes
            # 1. Calculate the derivative
            tmp1 .= -1im*(t_given[counter])*dHe*inbetween_states[counter, :]
            # 2. take the outer products and normalize
            xis[counter, :, :] = (adjoint(tmp1) .* inbetween_states[counter, : ] +
                               adjoint(inbetween_states[counter, : ]) .* tmp1) / dot(inbetween_states[counter, :], inbetween_states[counter, :])
            counter = counter + 1
            if counter > ntimes
                break
            end
        end
        return xis
    end
    # All the xis before the first jump can be handled like this:
    while (t_given[counter] < traj[counter_c].time) && (counter <= ntimes)
            # 1. Calculate the derivative
            tmp1 .= -1im*(t_given[counter])*dHe*inbetween_states[counter, :]
            # 2. normalize and add
            xis[counter, :, :] = (adjoint(tmp1) .* inbetween_states[counter, : ] +
                               adjoint(inbetween_states[counter, : ]) .* tmp1) / dot(inbetween_states[counter, :], inbetween_states[counter, :])
            counter = counter + 1
            if counter > ntimes
                break
            end
    end
    t_ = t_ + traj[counter_c].time
    counter_c = counter_c + 1
    # In between jumps
    while (counter_c <= njumps) && (counter <= ntimes)
        timeclick = traj[counter_c].time
        while (t_ < t_given[counter] < t_ + timeclick) && (counter <= ntimes)
             # 1. Calculate the derivative
             tmp1 .= -1im*(t_given[counter]-t_)*dHe*jump_states[counter_c, :] +
                     exp(-1im*(t_given[counter]-t_)*sys.Heff)*jump_dpsis[counter_c, :]
             xis[counter, :, :] = (adjoint(tmp1) .* inbetween_states[counter, : ] +
                               adjoint(inbetween_states[counter, : ]) .* tmp1) / dot(inbetween_states[counter, :], inbetween_states[counter, :])
             counter = counter + 1
             if counter > ntimes
                 break
             end
         end
       t_ = t_ + timeclick
       counter_c = counter_c + 1
    end
    # After all the jumps finished
    while counter <= ntimes
        tmp1 .= -1im*(t_given[counter]-t_)*dHe*inbetween_states[counter, : ] +
               exp(-1im*(t_given[counter]-t_)*sys.Heff)*jump_dpsis[end, :]
         xis[counter, :, :] = (adjoint(tmp1) .* inbetween_states[counter, : ] +
                               adjoint(inbetween_states[counter, : ]) .* tmp1) / dot(inbetween_states[counter, :], inbetween_states[counter, :])
        counter = counter + 1
    end
    return xis
end
