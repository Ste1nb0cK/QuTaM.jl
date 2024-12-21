sys = QuTaM.rf_sys
params = QuTaM.rf_params
data = run_trajectories(sys, params);
times = collect(LinRange(0,0.5,10));
data = QuTaM.sample_single_trajectory(sys, params);
states = QuTaM.states_at_jumps(data, sys, params.psi0)
@testset "States at Jumps" begin
    @test size(data) == size(states)
end
