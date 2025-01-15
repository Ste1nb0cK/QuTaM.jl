using QuTaM, LinearAlgebra, Test
@testset "Pure States" begin
    function identify_state(s::Vector{ComplexF64})
        s1 = [1.0+0im, 0]
        s2 = [0, 1.0+0im]
        if norm(s1-s) < 0.01
            return 1
        elseif norm(s2-s) < 0.01
            return 2
        else
            return 0
        end
    end
    sys = QuTaM.rdt_sys
    params = QuTaM.rdt_params
    traj = QuTaM.sample_single_trajectory(sys, params, 10)
    jump_dts = [click.time for click in traj]
    min_dt = extrema(jump_dts)[1]
    jump_times = cumsum(jump_dts)
    jump_labels = [click.label for click in traj]
    # Check that empty time vector returns empty vector
    t_given_empty = Vector{Float64}(undef, 0)
    states_empty = QuTaM.evaluate_at_t(t_given_empty, traj, sys, params.psi0)
    # Extremal case: empty t_given
    @test isempty(states_empty)
    # Check no-jump trajectory
    t_given_a = jump_times .+ min_dt/2
    # "No-jump trajectory"
    traj_empty = Vector{QuTaM.DetectionClick}(undef, 0)
    states_no_jumps = QuTaM.evaluate_at_t(t_given_a, traj_empty, sys, params.psi0)
    for k in 1:size(states_no_jumps)[1]
        @test identify_state(states_no_jumps[k, :]) == identify_state(params.psi0)
    end
    # Verify that evaluation just after the jumps coincides with the jumps
    states_between_jumps_a = QuTaM.evaluate_at_t(t_given_a, traj, sys, params.psi0)
    njumps = size(states_between_jumps_a)[1]
    # Evaluate just after jump and check unitarity
    for k in 1:njumps
        @test identify_state(states_between_jumps_a[k, :]) == jump_labels[k]
        @test abs(norm(states_between_jumps_a[k, :]) -1) < 0.01
    end
    # Evaluate just before jump
    t_given_b = jump_times .- min_dt/2
    states_between_jumps_b = QuTaM.evaluate_at_t(t_given_b, traj, sys, params.psi0)
    njumps = size(states_between_jumps_b)[1]
    for k in 1:njumps
        @test identify_state(states_between_jumps_b[k, :]) - 1 == (jump_labels[k] % 2)
        @test abs(norm(states_between_jumps_b[k, :]) -1) < 0.01
    end
end

@testset "Mixed States" begin
    function identify_state(s::Matrix{ComplexF64})
        excited = [[0, 0] [0, 1.0+0im]]
        ground = [[1.0+0im, 0] [0, 0]]
        if norm(ground-s) < 0.01
            return 1
        elseif norm(excited-s) < 0.01
            return 2
        else
            return 0
        end
    end
    rho0 = zeros(ComplexF64, 2, 2)
    # Initial condition: Completely mixed state
    rho0[2,2] = 0.5 # Initial condition
    rho0[1,1] = 0.5
    params = SimulParameters(rho0,
                                5.0, # Final time. Set very long so that all trajectories jump
                                10, # seed
                                1000, # Number of trajectories
                                10_000, # Number of samples in the finegrid
                                10.0, # Multiplier to use in the fine grid
                                1e-3 # Tolerance for passing Dark state test
                                )
    sys = QuTaM.rdt_sys
    traj = QuTaM.sample_single_trajectory(sys, params, 10)
    jump_dts = [click.time for click in traj]
    min_dt = extrema(jump_dts)[1]
    jump_times = cumsum(jump_dts)
    jump_labels = [click.label for click in traj]
    # Check that empty time vector returns empty vector
    t_given_empty = Vector{Float64}(undef, 0)
    states_empty = QuTaM.evaluate_at_t(t_given_empty, traj, sys, params.psi0)
    # Extremal case: empty t_given
    @test isempty(states_empty)
    # Check no-jump trajectory
    t_given_a = jump_times .+ min_dt/2
    # "No-jump trajectory"
    traj_empty = Vector{QuTaM.DetectionClick}(undef, 0)
    states_no_jumps = QuTaM.evaluate_at_t(t_given_a, traj_empty, sys, params.psi0;
                                          normalize=false)
    for k in 1:size(states_no_jumps)[1]
        @test norm(states_no_jumps[k, :, :] -
            exp(-1im*t_given_a[k]*sys.Heff) * params.psi0 * exp(1im*t_given_a[k]*adjoint(sys.Heff))) < 0.01
    end
    # Verify that evaluation just after the jumps is correct
    states_between_jumps_a = QuTaM.evaluate_at_t(t_given_a, traj, sys, params.psi0)
    njumps = size(states_between_jumps_a)[1]
    # Evaluate just after jump and check unitarity
    for k in 1:njumps
        @test identify_state(states_between_jumps_a[k, :, :]) == jump_labels[k]
        @test abs(tr(states_between_jumps_a[k, :, :]) -1) < 0.01
    end
end
