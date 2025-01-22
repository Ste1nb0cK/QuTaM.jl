# This test evaluates if the trace of the monitoring operator is working
using Test, QuTaM, LinearAlgebra
function analytical_contribution(t, traj::Trajectory, gamma::Float64, n)
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
    #print(m)
    return m/gamma - sum(jump_times[1:m]) -(n-m)*t

end

qubit_sys = QuTaM.rd_sys
qubit_params = SimulParameters(QuTaM.rd_psi0,
    3.0, # Final time. Set very long so that all trajectories jump
    1, # seed
    100, # Number of trajectories
    50_000, # Number of samples in the finegrid
    10.5, # Multiplier to use in the fine grid
    1e-3 # Tolerance for passing Dark state test
)
H_parametrized_qubit = (delta::Float64, gamma::Float64) -> (0.5*delta*QuTaM.sigma_z)::Matrix{ComplexF64}
L_parametrized_qubit = (delta::Float64, gamma::Float64) -> (sqrt(gamma)*QuTaM.sigma_m)::Matrix{ComplexF64}
Heff_parametrized_qubit = QuTaM.GetHeffParametrized(H_parametrized_qubit, [L_parametrized_qubit])

qubit_trajectories = run_trajectories(qubit_sys, qubit_params);

ntimes = 100
t_given = collect(LinRange(0, qubit_params.tf, ntimes));

qubit_xi_sample = Array{ComplexF64}(undef, qubit_sys.NLEVELS, qubit_sys.NLEVELS, ntimes, qubit_params.ntraj)
for n in 1:qubit_params.ntraj
    qubit_xi_sample[:, :, :, n] = MonitoringOperator(t_given, qubit_sys, Heff_parametrized_qubit,
                                                    [L_parametrized_qubit], qubit_trajectories[n], qubit_params.psi0,
                           [QuTaM.rd_deltaomega, QuTaM.rd_gamma], [0.0, QuTaM.rd_gamma/100])
end

contribution_qubit_sample = Array{Float64}(undef, ntimes, qubit_params.ntraj)
for n in 1:qubit_params.ntraj
    for k in 1:ntimes
       contribution_qubit_sample[k, n] = real(tr(qubit_xi_sample[:, :, k, n]))
    end
end
qubit_EPS = 0.1 # Difference tolerance for the test
#f_analytical_qubit(t) = analytical_contribution(t, trajectories[10], qutam.rd_gamma, sys.nlevels-1)
##################### Single jump Trajectories #######
@testset "Single Jump Trajectories" begin
    for k in 1:qubit_params.ntraj
        for t in 1:ntimes
            flag = abs(contribution_qubit_sample[t, k] - analytical_contribution(t_given[t], qubit_trajectories[k], QuTaM.rd_gamma, qubit_sys.NLEVELS-1)) < qubit_EPS
            flag || @warn "Failure at trajectory=$k, t=$t"
            @test flag
        end
    end
end
