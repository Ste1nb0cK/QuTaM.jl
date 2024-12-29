using QuTaM, LinearAlgebra
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
@testset "Empty t_given" begin
   @test isempty(states_empty)
end
# Check no-jump trajectory
t_given_a = jump_times .+ min_dt/2
@testset "No-jump trajectory" begin
    traj_empty = Vector{QuTaM.DetectionClick}(undef, 0)
    states_no_jumps = QuTaM.evaluate_at_t(t_given_a, traj_empty, sys, params.psi0)
    for k in 1:size(states_no_jumps)[1]
        @test identify_state(states_no_jumps[k, :]) == identify_state(params.psi0)
    end
end
# Verify that evaluation just after the jumps coincides with the jumps
@testset "States after jumps " begin
    states_between_jumps_a = QuTaM.evaluate_at_t(t_given_a, traj, sys, params.psi0)
    njumps = size(states_between_jumps_a)[1]
   # Evaluate just after jump and check unitarity
    for k in 1:njumps
        @test identify_state(states_between_jumps_a[k, :]) == jump_labels[k]
        @test abs(norm(states_between_jumps_a[k, :]) -1) < 0.01
    end
end
# Evaluate just before jump
@testset "States before jumps " begin
    t_given_b = jump_times .- min_dt/2
    states_between_jumps_b = QuTaM.evaluate_at_t(t_given_b, traj, sys, params.psi0)
    njumps = size(states_between_jumps_b)[1]
    for k in 1:njumps
        @test identify_state(states_between_jumps_b[k, :]) - 1 == (jump_labels[k] % 2)
        @test abs(norm(states_between_jumps_b[k, :]) -1) < 0.01
    end
end
