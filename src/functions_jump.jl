# Function that returns a view of the given array, with the last index fixed at k
function fixlastindex(array::Array{ComplexF64}, k::Int64)
    indices = ntuple(d -> d == ndims(array) ? k : Colon(), ndims(array))
    # Add a singleton of dimension 1
    return view(array, indices...)
end

function calculatewtdweights!(W::Array{Float64}, Qs::Array{ComplexF64}, psi::Vector{ComplexF64}, params::SimulParameters)
    @inbounds @simd for k in 1:params.nsamples
           W[k] = real(dot(psi, Qs[:, :, k], psi)) # dot product without storing A*x. THIS IS THE KEY FOR SPEED
        end
end

function calculatewtdweights!(W::Array{Float64}, Qs::Array{ComplexF64}, psi::Matrix{ComplexF64}, params::SimulParameters)
    @inbounds @simd for k in 1:params.nsamples
           W[k] = real(tr(Qs[:, :, k] * psi))
        end
end

function calculatechannelweights!(P::Vector{Float64}, psi::Vector{ComplexF64}, sys::System)
    aux_P = real(dot(psi, sys.J * psi))
    @inbounds @simd for k in 1:sys.NCHANNELS
        P[k] = norm(sys.Ls[k]*psi)^2
    end
    P .= P / aux_P
end

function calculatechannelweights!(P::Vector{Float64}, psi::Matrix{ComplexF64}, sys::System)
    aux_P = real(dot(psi, sys.J * psi))
    @inbounds @simd for k in 1:sys.NCHANNELS
        P[k] = real(tr(sys.LLs[k]*psi))
    end
    P .= P / aux_P
end

# Pure state versions
function prejumpupdate!(V::Matrix{ComplexF64}, psi::Vector{ComplexF64}; normalize=false)
    psi .= V * psi
    if normalize
        psi .= psi/ norm(psi)
    end
end

function prejumpupdate!(psi::Vector{ComplexF64}, V::Matrix{ComplexF64},
                        psi0::Union{Vector{ComplexF64}, SubArray{ComplexF64}}; normalize=false)
    psi .= V * psi0
    if normalize
        psi .= psi/ norm(psi)
    end
end

# Mixed state versions

function prejumpupdate!(psi::Matrix{ComplexF64}, V::Matrix{ComplexF64},
                        psi0::Union{Matrix{ComplexF64}, SubArray{ComplexF64}}; normalize=false)
    psi .= V * psi0 * adjoint(V)
    if normalize
        psi .= psi/ tr(psi)
    end
end

function prejumpupdate!(V::Matrix{ComplexF64}, psi::Matrix{ComplexF64}; normalize=false)
    psi .= V * psi * adjoint(V)
    if normalize
        psi .= psi/ tr(psi)
    end
end

function postjumpupdate!(L::Matrix{ComplexF64}, psi::Vector{ComplexF64}; normalize=true)
        psi .= L*psi # State without normalization
        if normalize
            psi .= psi / norm(psi)
        end
end

function postjumpupdate!(L::Matrix{ComplexF64}, psi::Matrix{ComplexF64}; normalize=true)
        psi .= L*psi*adjoint(L) # State without normalization
        if normalize
            psi .= psi / tr(psi)
        end
end


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
function run_singletrajectory(sys::System, params::SimulParameters,
    W::Vector{Float64}, P::Vector{Float64}, ts::Vector{Float64},
    Qs::Array{ComplexF64}, Vs::Array{ComplexF64}; seed::Int64 = 1, isrenewal=false)
    Random.seed!(seed)
    channel = 0
    traj = Vector{DetectionClick}()
    psi = copy(params.psi0)
    t::Float64 = 0
    channel = 0
    # Run the trajectory
    calculatewtdweights!(W, Qs, psi, params)
    while t < params.tf
        #Sample jump time and  move state to pre-jump state
        tau_index = StatsBase.sample(1:params.nsamples, StatsBase.weights(W))
        t = ts[tau_index] + t
        prejumpupdate!(Vs[:, :, tau_index], psi)
        # Sample jump channel
        calculatechannelweights!(P, psi, sys)
        channel = StatsBase.sample(1:sys.NCHANNELS, StatsBase.weights(P))
        # State update
        postjumpupdate!(sys.Ls[channel], psi)
        push!(traj, DetectionClick(ts[tau_index], channel))
        # Sample WTD
        if !isrenewal
            calculatewtdweights!(W, Qs, psi, params)
            if sum(W) < params.eps
                break
            end
        end
    end
    return traj
end

function writestate!(states::Array{ComplexF64},
                     psi::Union{Vector{ComplexF64}, Matrix{ComplexF64}}, counter::Int64)
             fixlastindex(states, counter) .= psi
end



# function writestate!(states::Array{ComplexF64}, psi::Vector{ComplexF64}, counter::Int64)
            # states[:, counter] .= psi
# end

# function writestate!(states::Array{ComplexF64}, psi::Matrix{ComplexF64}, counter::Int64)
            # states[:, :, counter] .= psi
# end


############ Evaluation at given times #######################
"""
    statesat_jumps(traj::Trajectory, sys::System,
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
function states_atjumps(traj::Trajectory, sys::System,
                      psi0::Union{Vector{ComplexF64}, Matrix{ComplexF64}}; normalize::Bool=true)
    njumps = size(traj)[1]
    if isa(psi0, Vector{ComplexF64})
        states = Array{ComplexF64}(undef, sys.NLEVELS,  njumps)
    elseif isa(psi0, Matrix{ComplexF64})
        states = Array{ComplexF64}(undef, sys.NLEVELS, sys.NLEVELS, njumps)
    end
    psi = copy(psi0)
    jump_counter = 1
    for click in traj
        prejumpupdate!(exp(-1im*(click.time)*sys.Heff), psi)
        postjumpupdate!(sys.Ls[click.label], psi; normalize=normalize)
        writestate!(states, psi, jump_counter)
        jump_counter = jump_counter + 1
    end
    return states
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

function states_att(t_given::Vector{Float64}, traj::Trajectory, sys::System,
                       psi0::Union{Vector{ComplexF64}, Matrix{ComplexF64}};
                       normalize::Bool=true)
    # Special case: if the time array is empty, return an empty array
    if isempty(t_given)
        return Array{ComplexF64}(undef, 0, 0) # empty 2 dimensional array
    end
    psi = copy(psi0)
    ntimes = size(t_given)[1]
    jump_states = states_atjumps(traj, sys, psi0; normalize=normalize)
    njumps = size(traj)[1]
    t_ = 0
    counter = 1
    counter_c = 1
    # states = Array{ComplexF64}(undef, sys.NLEVELS, ntimes)
    if isa(psi0, Vector{ComplexF64})
        states = Array{ComplexF64}(undef, sys.NLEVELS,  ntimes)
    elseif isa(psi0, Matrix{ComplexF64})
        states = Array{ComplexF64}(undef, sys.NLEVELS, sys.NLEVELS, ntimes)
    end
    # Edge case: if the trajectory is empty, evaluate exponentials and return
    if isempty(traj)
        while counter <= ntimes
            prejumpupdate!(psi, exp(-1im*(t_given[counter])*sys.Heff), psi0;
                           normalize=normalize)
            # fixlastindex(states, counter)
            writestate!(states, psi, counter)
            counter = counter + 1
            if counter > ntimes
                break
            end
        end
        return states
    end
    # All the states before the first jump can be handled like this:
    while (t_given[counter] < traj[counter_c].time) && (counter <= ntimes)
            prejumpupdate!(psi, exp(-1im*(t_given[counter])*sys.Heff), psi0;
                           normalize=normalize)
            writestate!(states, psi, counter)
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
             prejumpupdate!(psi, exp(-1im*(t_given[counter] - t_)*sys.Heff),
                            fixlastindex(jump_states, counter_c-1); normalize=normalize)
             writestate!(states, psi, counter)
             counter = counter + 1
             if counter > ntimes
                 break
             end
         end
       t_ = t_ + timeclick
       counter_c = counter_c + 1
    end

    while counter <= ntimes
        prejumpupdate!(psi, exp(-1im*(t_given[counter] - t_)*sys.Heff),
                       fixlastindex(jump_states, njumps); normalize=normalize)
        writestate!(states, psi, counter)
        counter = counter + 1
    end
    return states
end


