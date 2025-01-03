################ SETUP FOR RESONANCE FLUORESENCE ##################
# Same thing as Radiative Damping, but now there Hamiltonian is different
rf_EPS = 1e-5 # Tolerance for the distance respect to the Frobenious norm
rf_delta= 1.43/2 # Detuning
rf_omega = 1.3/2 # Rabi Frequency
rf_gamma = 0.55
rf_H = rf_delta*sigma_z + rf_omega*sigma_x

rf_psi0 = zeros(ComplexF64, 2)
rf_psi0[1] = 1 # Initial condition

rf_sys = System(rf_H, # Hamiltonian
[sqrt(rf_gamma)*sigma_m]) #Jump Operators
rf_params = SimulParameters(rf_psi0,
    400.0, # Final time. Set very long so that all trajectories jump
    1, # seed
    250, # Number of trajectories
    75_000, # Number of samples in the finegrid
    4.0, # Multiplier to use in the fine grid
    1e-3 # Tolerance for passing Dark state test
)

# Distribution class for
