# This test evaluates if the trace of the monitoring operator is working
using Test, BackAction, LinearAlgebra
####################### CASE: HAMILTONIAN COMMUTES WITH ITS DERIVATIVE
# Consider a radiative damping for the case of an N-level system.
# The logarithmic derivatives of the trajectory distribution
@testset "Case: Hamiltonian commutes with its derivative" begin
    function analytical_contribution(t, traj::Trajectory, gamma::Float64, n)
        # Edge case: no click trajectory
        if isempty(traj)
            return -n*t
        end
        jump_times = cumsum([click.time for click in traj])
        m = 0 # number of jumps
        # Find which was the last jump that occured
        for t_jump in jump_times
            if t_jump > t
                break
            end
            m = m +1
        end
        if m == 0
            return -n*t
        end
        return m/gamma - sum(jump_times[1:m]) -(n-m)*t

    end

    ############# System Definition
    NLEVELS = 5
    # Construct the destruction operator
    a = zeros(ComplexF64, NLEVELS, NLEVELS)
    for k in 1:NLEVELS
        for m in 1:NLEVELS
            if k == m + 1
                a[m, k] = sqrt(k-1)
            end
        end
    end

    EPS = 1e-5 # Tolerance
    GAMMA = BackAction.rd_gamma

    H = 0.5*BackAction.rd_deltaomega*adjoint(a)*a
    L = sqrt(GAMMA)*a
    sys = System(H, [L]);

    PSI0 = zeros(ComplexF64, NLEVELS)
    PSI0[end] = 1 # Initial condition

    params = SimulParameters(PSI0,
                             3.0, # Final time. Set very long so that all trajectories jump
                             2, # seed
                             100, # Number of trajectories
                             50_000, # Number of samples in the finegrid
                             10.5, # Multiplier to use in the fine grid
                             1e-3 # Tolerance for passing Dark state test
                             );

    # Parametrization
    H_parametrized = (delta::Float64, gamma::Float64) -> (0.5*delta*adjoint(a)*a)::Matrix{ComplexF64}
    L_parametrized = (delta::Float64, gamma::Float64) -> (sqrt(gamma)*a)::Matrix{ComplexF64}
    He_parametrized = BackAction.getheff_parametrized(H_parametrized, [L_parametrized])
    trajectories = run_trajectories(sys, params)

    ntimes = 1000
    t_given = collect(LinRange(0, params.tf, ntimes));

    EPS = 0.1 # Difference tolerance for the test

    xi_sample = Array{ComplexF64}(undef, sys.NLEVELS, sys.NLEVELS, ntimes, params.ntraj)
    theta = [BackAction.rd_deltaomega, BackAction.rd_gamma]
    dtheta = [0.0, BackAction.rd_gamma/100]
    for n in 1:params.ntraj
        xi_sample[:, :, :, n] = monitoringoperator(t_given, sys, He_parametrized, [L_parametrized], trajectories[n], params.psi0,
                                                   theta, dtheta)
    end
    contribution_sample = Array{Float64}(undef, ntimes, params.ntraj)
    for k in 1:params.ntraj
        for tn in 1:ntimes
            contribution_sample[tn, k] = real(tr(xi_sample[:, :, tn, k]))
        end
    end
    # Check that the contribution from the monitoring operator coincides with the analytical one
    for k in 1:params.ntraj
        for t in 1:ntimes
            flag = abs(
                contribution_sample[t, k] - analytical_contribution(t_given[t],
                                                                    trajectories[k],
                                                                    BackAction.rd_gamma,
                                                                    sys.NLEVELS-1)) < EPS
            flag || @warn "Failure at trajectory=$k, t=$t"
            @test flag
        end
    end
end
