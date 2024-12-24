################ SETUP FOR RESONANCE FLUORESENCE ##################
# Same thing as Radiative Damping, but now there Hamiltonian is different
rf_EPS = 1e-5 # Tolerance for the distance respect to the Frobenious norm
rf_delta= 1.0 # Detuning
rf_omega = 1.1 # Rabi Frequency
rf_gamma = 10.0
rf_H = 0.5*rf_delta*sigma_z + 0.5*rf_omega*sigma_x
#WARNING: THE DERIVED OPERATORS ARE BROKEN
# rf_L = sqrt(rf_gamma) * sigma_m
# rf_J = rf_gamma * [[0,0] [0,1.0+0im]]
# rf_He = [[rf_deltaomega/2, 0.0] [0.0, 0.5*(-rf_deltaomega - 1im*rf_gamma) ]]

psi0 = zeros(ComplexF64, 2)
psi0[2] = 1 # Initial condition

rf_sys = System(rf_H, # Hamiltonian
[sqrt(rf_gamma)*sigma_m]) #Jump Operators
rf_params = SimulParameters(psi0,
    10.0, # Final time. Set very long so that all trajectories jump
    1, # seed
    1000, # Number of trajectories
    5_000, # Number of samples in the finegrid
    5.0, # Multiplier to use in the fine grid
    1e-3 # Tolerance for passing Dark state test
)
