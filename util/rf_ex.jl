############## SETUP FOR RESONANCE FLOURESCENE
rf_gamma = 1.5
rf_n = 1
rf_gamma1 = (rf_n+1)*rf_gamma
rf_gamma2 = (rf_n)*rf_gamma
rf_sys = System( zeros(ComplexF64, 2, 2), # Hamiltonian
    [sqrt(rf_gamma1)*QuTaM.sigma_m, sqrt(rf_gamma2)*QuTaM.sigma_p ]) #Jump Operators
#### 2. Create the simulation parameters instance
rf_psi0 = zeros(ComplexF64, 2)
rf_psi0[2] = 1 # Initial condition
rf_params= SimulParameters(psi0,
    5.0, # Final time. Set very long so that all trajectories jump
    10, # seed
    1000, # Number of trajectories
    10_000, # Number of samples in the finegrid
    10, # Multiplier to use in the fine grid
    1e-3 # Tolerance for passing Dark state test
)
