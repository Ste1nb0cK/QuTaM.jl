############# Precomputation routine ######################
# Input:
# 1. From the System: J and Heff
# 2. From the simulation parameters: tf, nsamples and the multiplier
# Output: (ts, Qs), a tuple with the finegrid and the precomputed values of the transformed J
function precompute(J::Matrix{ComplexF64}, Heff::Matrix{ComplexF64},
        tf::Float64, nsamples::Int64, multiplier::Float64)
        ts = LinRange(0, multiplier*tf, nsamples)
        Qs = Vector{Matrix{ComplexF64}}(undef, nsamples)
        for k in 1:nsamples
            expm = exp(-1im*ts[k]*Heff)
            Qs[k] = expm * J * adjoint(expm)
        end
        return ts, Qs
end
############# Single Trajectory Routine ######################
# Input:
# 1. From the System: Ls and Heff
# 2. From the Simulation Parameters: psi0, tf, nsamples, dt and eps
# Output:
# A trajectory object with the data of the trajectory
function run_single_trajectory(sys::System, params::SimulParameters)
    # System
    J = sys.J
    Heff = sys.Heff
    Ls = sys.Ls
    # Simulation parameters
    tf = params.tf
    nsamples = params.nsamples
    dt = params.dt
    eps = params.eps
    # Random number generator
    Random.seed!(params.seed)
    # Store
    W = zeros(Float64, nsamples) # Vector to store the weights of the fine grid
    times = Vector{Float64}()
    labels = Vector{Int64}()
    states = Vector{Vector{ComplexF64}}()

   # 1. Set Initial condition
    psi = params.psi0
    t = 0
    # 2. Precomputing
    ts, Qs =  precompute(J, Heff, tf, nsamples, params.multiplier)
    # 3. Run the trajectory
    while t < tf
    # 1. Calculate the WTD for the state, these act as weights
        for k in 1:nsamples
           W[k] = real(dot(conj.(psi), Qs[k]*psi))
        end
        # 1.a We must verify if we got a dark state, that can be cheked by
        # looking at the normalization of the QTD
        if abs(sum(W)*dt - 1) > eps
            break
        end
        # 2. Sample jump time
        tau = StatsBase.sample(ts, StatsBase.weights(W))
        t = tau + t
        if t>tf # If the next jump happens after tf, stop
            break
        end
        psi_tilde = exp(-1im*tau*Heff) * Ls[1]*psi # State without normalization
        psi = psi_tilde / norm(psi_tilde)
        push!(states, psi)
        push!(labels, 1)
        push!(times, t)
    end
    return Trajectory(times, states, labels)
end
############## Multitrajectory Routine ###########################
function run_trajectories(sys::System, params::SimulParameters)
    # System
    J = sys.J
    Heff = sys.Heff
    Ls = sys.Ls
    # Simulation parameters
    tf = params.tf
    nsamples = params.nsamples
    dt = params.dt
    eps = params.eps
    ntraj = params.ntraj
   # Store
    W = zeros(Float64, nsamples) # Vector to store the weights of the fine grid
    times = Vector{Float64}()
    labels = Vector{Int64}()
    states = Vector{Vector{ComplexF64}}()
    data = Vector{Trajectory}(undef, ntraj)
    # Precomputing
    ts, Qs =  precompute(J, Heff, tf, nsamples, params.multiplier)
    for k in 1:ntraj
        # Set Initial condition
        psi = params.psi0
        t = 0
        # Random number generator
        Random.seed!(params.seed + k)
       #  Run the trajectory
        while t < tf
        # 1. Calculate the WTD for the state, these act as weights
            for k in 1:nsamples
               W[k] = real(dot(conj.(psi), Qs[k]*psi))
            end
            # 1.a We must verify if we got a dark state, that can be cheked by
            # looking at the normalization of the QTD
            if abs(sum(W)*dt - 1) > eps
                break
            end
            # 2. Sample jump time
            tau = StatsBase.sample(ts, StatsBase.weights(W))
            t = tau + t
            if t>tf # If the next jump happens after tf, stop
                break
            end
            psi_tilde = exp(-1im*tau*Heff) * Ls[1]*psi # State without normalization
            psi = psi_tilde / norm(psi_tilde)
            push!(states, psi)
            push!(labels, 1)
            push!(times, t)
        end
            data[k] =  Trajectory(copy(times), copy(states), copy(labels))
            empty!(states)
            empty!(labels)
            empty!(times)
    end
    return data
end
