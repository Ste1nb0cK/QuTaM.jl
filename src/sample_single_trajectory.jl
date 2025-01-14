
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
    psi = Matrix{ComplexF64}(undef, sys.NLEVELS, sys.NLEVELS)
    W = Vector{Float64}(undef, params.nsamples)
    P = Vector{Float64}(undef, sys.NCHANNELS)
    data = run_single_trajectory(sys, params,
                                W, P, psi, ts, Qs, seed = seed)
    return data
end
