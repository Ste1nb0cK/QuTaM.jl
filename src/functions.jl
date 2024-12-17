############# Precomputation routine ######################
# Input:
# 1. From the System: J and Heff
# 2. From the simulation parameters:  nsamples
# Output: None
# Action: Modifies Qs to store the precomputed Q(ts)
function precompute!(sys::System, nsamples::Int64,
         ts::Vector{Float64}, Qs::Vector{Matrix{ComplexF64}})
        for k in 1:nsamples
            expm = exp(-1im*ts[k]*sys.Heff)
            Qs[k] = expm * sys.J * adjoint(expm)
        end
    return
end
############# Single Trajectory Routine ######################
# Input:
# 1. From the System: Ls and Heff
# 2. From the Simulation Parameters it uses: psi0, tf, nsamples, dt, eps and seed
# 3. Extra:  vectors W to store the weights of the jump times, a vector
# P to store the weights on the chhanels, (ts, Qs) from the precomputing and the seed
# of the single trajectory BEWARE: IT MIGHT NOT COINCIDE WITH THAT OF THE SIMULATION PARAMETERS,
# a vector psi to store the current state
# Output:
# A trajectory object with the data of the trajectory
function run_single_trajectory(
    sys::System,
    params::SimulParameters,
    W::Vector{Float64}, P::Vector{Float64}, psi::Vector{ComplexF64}, ts::Vector{Float64},
    Qs::Vector{Matrix{ComplexF64}}; seed::Int64 = 1)
    # Random number generator
    Random.seed!(seed)
    # Store
    times = Vector{Float64}()
    labels = Vector{Int64}()
    states = Vector{Vector{ComplexF64}}()

   # 1. Set Initial condition
    psi .= params.psi0
    t::Float64 = 0
    channel = 0
    # 2. Run the trajectory
    while t < params.tf
        # If the probability of no jump is above the tolerance, declare dark state
        q0 = norm(exp(-1im*params.multiplier*params.tf*sys.Heff)*psi)
        if q0^2 > params.eps
            break
        end
        # Calculate the WTD for the state, these act as weights
        for k in 1:params.nsamples
           W[k] = real(dot(psi, Qs[k]*psi))
        end
        # 2. Sample jump time
        tau = StatsBase.sample(ts, StatsBase.weights(W))
        t = tau + t
        if t > params.tf # If the next jump happens after tf, stop
            break
        end
        psi .= exp(-1im*tau*sys.Heff) * psi
        # 3. Sample the channel
        aux_P = real(dot(psi, sys.J * psi))
        for k in 1:sys.NCHANNELS
            P[k] = norm(sys.Ls[k]*psi)^2
        end
        P .= P / aux_P
        channel = StatsBase.sample(1:sys.NCHANNELS, StatsBase.weights(P))
        psi .= sys.Ls[channel]*psi # State without normalization
        psi .= psi / norm(psi)
        push!(states, psi)
        push!(labels, channel)
        push!(times, t)
    end
    return Trajectory(times, states, labels)
end
############## Multitrajectory Routine ###########################
function run_trajectories(sys::System, params::SimulParameters)
    ## Precomputing
    ts = collect(LinRange(0, params.multiplier*params.tf, params.nsamples))
    Qs = Vector{Matrix{ComplexF64}}(undef, params.nsamples)
    precompute!(sys, params.nsamples, ts, Qs)
    # Running the trajectory
    psi = Vector{ComplexF64}(undef, sys.NLEVELS)
    W = Vector{Float64}(undef, params.nsamples)
    P = Vector{Float64}(undef, sys.NCHANNELS)
    data = Vector{Trajectory}(undef, params.ntraj)
    for k in 1:params.ntraj
        data[k] = run_single_trajectory(sys, params,
                                        W, P, psi, ts, Qs, seed = params.seed + k)
    end
    return data
end
