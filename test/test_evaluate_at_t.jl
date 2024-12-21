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
@testset "Inbetween Jump Evaluation" begin
    sys = QuTaM.rdt_sys
    params = QuTaM.rdt_params
    traj = QuTaM.sample_single_trajectory(sys, params, 10)
    jump_dts = [click.time for click in traj]
    min_dt = extrema(jump_dts)[1]
    jump_times = cumsum(jump_dts)
    jump_labels = [click.label for click in traj]
    # Verify that evaluation just after the jumps coincides with the jumps
    t_given = jump_times .+ min_dt/2
    states_after_jumps = QuTaM.states_at_jumps(traj, sys, params.psi0)
    njumps = size(states_after_jumps)[1]
    # Evaluate just after jump
    for k in 1:njumps
    @test identify_state(states_after_jumps[k]) == jump_labels[k]
    end
end
