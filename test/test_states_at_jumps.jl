using QuTaM, LinearAlgebra, Test

@testset "Pure States" begin
    ############### Check sizes match for multiple jumps
    sys = QuTaM.rdt_sys
    params = QuTaM.rdt_params
    for seed in 1:3
       traj = QuTaM.sample_single_trajectory(sys, params, seed);
       states = QuTaM.states_at_jumps(traj, sys, params.psi0)
       @test size(traj)[1] == size(states)[2]
   end
    ################### Check the correct states are obtained
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
    counter = 1
    traj = QuTaM.run_trajectories(sys, params)[1]
    states = QuTaM.statesatjumps(traj, sys, params.psi0)
    for click in traj
        @test identify_state(states[:, counter]) == click.label
        counter = counter + 1
    end
    ############## Edge case: no jumps
    psi0 = zeros(ComplexF64, 2)
    psi0[2] = 1 # Initial condition
    sys = System(QuTaM.rd_H, # Hamiltonian
                 [sqrt(QuTaM.rd_gamma)*QuTaM.sigma_m]) #Jump Operators
    params = SimulParameters(psi0,
                             1e-4, # Final time. Set very long so that all trajectories jump
                             1, # seed
                             20, # Number of trajectories
                             3_000, # Number of samples in the finegrid
                             3.0, # Multiplier to use in the fine grid
                             1e-3 # Tolerance for passing Dark state test
                             )
    traj = Vector{DetectionClick}(undef, 0) #QuTaM.sample_single_trajectory(sys, params, params.seed)
    states = QuTaM.states_at_jumps(traj, sys, params.psi0)
   @test isempty(states)
end


@testset "Mixed States" begin
    ############### Check sizes match for multiple jumps
    sys = QuTaM.rdt_sys
    #### Same thing as params in rd_temperature. The difference is the initial state
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

    for seed in 1:3
       traj = QuTaM.sample_single_trajectory(sys, params, seed);
       states = QuTaM.statesatjumps(traj, sys, params.psi0)
       @test size(traj)[1] == size(states)[3]
   end
    ################### Check the correct states are obtained
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
    counter = 1
    traj = QuTaM.sample_single_trajectory(sys, params, params.seed)
    states = QuTaM.statesatjumps(traj, sys, params.psi0)
    for click in traj
        @test identify_state(states[:, :, counter ]) == click.label
        counter = counter + 1
    end
    ############## Edge case: no jumps
    rho0 = zeros(ComplexF64, 2, 2)
    rho0[2, 2] = 0.5 # Initial condition
    rho0[1, 1] = 0.5 # Initial condition
    sys = System(QuTaM.rd_H, # Hamiltonian
                 [sqrt(QuTaM.rd_gamma)*QuTaM.sigma_m]) #Jump Operators
    params = SimulParameters(rho0,
                             1e-4, # Final time. Set very long so that all trajectories jump
                             1, # seed
                             20, # Number of trajectories
                             3_000, # Number of samples in the finegrid
                             3.0, # Multiplier to use in the fine grid
                             1e-3 # Tolerance for passing Dark state test
                             )
    traj = Vector{DetectionClick}(undef, 0) #QuTaM.sample_single_trajectory(sys, params, params.seed)
    states = QuTaM.states_at_jumps(traj, sys, params.psi0)
   @test isempty(states)
end
