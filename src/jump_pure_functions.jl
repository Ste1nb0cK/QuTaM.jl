############# Single Trajectory Routine ######################
"""
    run_single_trajectories(sys::System, params::SimulParameters,
     W::Vector{Float64}, P::Vector{Float64}, psi::Vector{ComplexF64},
     ts::Vector{Float64}, Qs::Vector{Matrix{ComplexF64}};
      seed::Int64 = 1) -> Trajectory

Sample a single trajectory from the system and parameters using the Gillipsie algorithm
. This is inteded to be used by `run_trajectories`

# Requiered Arguments:
- `sys::System`: System of interest
- `params::SimulParameters`: simulation parameters
- `W::Vector{Float64}`: to store the weights over the fine grid.
- `P::Vector{Float64}`: to store the weights over the channels
- `psi::Vector{ComplexF64}`: to store the current state vector
- `ts::Vector{Float64}`: the finegrid of waiting times
- `Qs::Vector{Matrix{ComplexF64}}`: 3-dimensional array for storing, dimensions must be (sys.NLEVELS, sys.NLEVELS, nsamples)

# Optional Arguments:
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
    W::Vector{Float64}, P::Vector{Float64}, ts::Vector{Float64},
    Qs::Array{ComplexF64}, Vs::Array{ComplexF64}; seed::Int64 = 1)
    # Random number generator
    Random.seed!(seed)
    traj = Vector{DetectionClick}()
    psi = copy(params.psi0)
    # temp_psi = similar(psi)
    t::Float64 = 0
    channel = 0
    # Run the trajectory
    while t < params.tf
        # Calculate the probability at infinity
        @inbounds @simd for k in 1:params.nsamples
           W[k] = real(dot(psi, Qs[:, :, k], psi)) # dot product without storing A*x. THIS IS THE KEY FOR SPEED
        end
        if sum(W) < params.eps
            break
        end
        # 2. Sample jump time
        # tau = StatsBase.sample(ts, StatsBase.weights(W))
        tau_index = StatsBase.sample(1:params.nsamples, StatsBase.weights(W))
        t = ts[tau_index] + t
        # psi .= exp(-1im*tau*sys.Heff) * psi
        psi .= Vs[:, :, tau_index] * psi
        # 3. Sample the channel
        aux_P = real(dot(psi, sys.J * psi))
        @inbounds @simd for k in 1:sys.NCHANNELS
            P[k] = norm(sys.Ls[k]*psi)^2
        end
        P .= P / aux_P
        channel::Int64 = StatsBase.sample(1:sys.NCHANNELS, StatsBase.weights(P))
        psi .= sys.Ls[channel]*psi # State without normalization
        psi .= psi / norm(psi)
        push!(traj, DetectionClick( ts[tau_index], channel))
    end
    return traj
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

The dimensions of the returned array `s` are `(sys.NLEVELS, size(traj))`,
so to recover the state vector at the ``n``-th jump one would do `s[:, n]`.
"""
function states_at_jumps(traj::Trajectory, sys::System,
                      psi0::Vector{ComplexF64}; normalize::Bool=true)
    njumps = size(traj)[1]
    # states = Vector{Vector{ComplexF64}}(undef, njumps)
    states = Array{ComplexF64}(undef, sys.NLEVELS,  njumps)
    psi = copy(psi0)
    jump_counter = 1
    if normalize
        for click in traj
            psi .= sys.Ls[click.label] * exp(-1im*(click.time)*sys.Heff) * psi
            psi .= psi/norm(psi)
            states[:, jump_counter] = psi[:]
            jump_counter = jump_counter + 1
        end
        return states
    else
        for click in traj
            psi .= sys.Ls[click.label] * exp(-1im*(click.time)*sys.Heff) * psi
            states[:, jump_counter] = psi[:]
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

# Required Arguments
- `t_given::Vector{Float64}`: times at which the trajectory is to be evalauted
- `traj::Trajectory`: the trajectory
- `sys::System`: the system to which the trajectory corresponds
- `psi0::Vector{ComplexF64}`: the initial state of the trajectory
# Optional Arguments
- `normalize::Bool`: whether to normalize the states or not, true by default
# Returns
A complex two-dimensional array whose rows contain the states.
"""

function evaluate_at_t(t_given::Vector{Float64}, traj::Trajectory, sys::System,
                       psi0::Vector{ComplexF64};
                       normalize::Bool=true)
    # Special case: if the time array is empty, return an empty array
    if isempty(t_given)
        return Array{ComplexF64}(undef, 0, 0) # empty 2 dimensional array
    end
    psi = copy(psi0)
    ntimes = size(t_given)[1]
    jump_states = states_at_jumps(traj, sys, psi0; normalize=normalize)
    njumps = size(traj)[1]
    t_ = 0
    counter = 1
    counter_c = 1
    states = Array{ComplexF64}(undef, sys.NLEVELS, ntimes)

    # Edge case: if the trajectory is empty, evaluate exponentials and return
    if isempty(traj)
        while counter <= ntimes
            psi .= exp(-1im*(t_given[counter])*sys.Heff) * psi0
            if normalize
                psi .= psi/norm(psi)
            end
            states[:, counter] = psi[:]
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
            states[:, counter] = psi[:]
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
             psi .= exp(-1im*(t_given[counter] - t_)*sys.Heff) * jump_states[:, counter_c-1]
             if normalize
             psi .= psi/norm(psi)
             end
             states[:, counter] .= psi[:]
             counter = counter + 1
             if counter > ntimes
                 break
             end
         end
       t_ = t_ + timeclick
       counter_c = counter_c + 1
    end

    while counter <= ntimes
        psi .= exp(-1im*(t_given[counter] - t_)*sys.Heff) * jump_states[:, njumps]
        if normalize
            psi .= psi/norm(psi)
        end
        states[ :, counter] = psi[:]
        counter = counter + 1
    end
    return states
end


