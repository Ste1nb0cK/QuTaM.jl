EPS = 1e-5 # Tolerance for the distance respect to the Frobenious norm
deltaomega = 10.0
gamma = 3.0
sigma_z = [[1.0+0im, 0] [0, -1]]
sigma_m = [[0.0+0im, 0] [1, 0]]
H = 0.5*deltaomega * sigma_z
L = sqrt(gamma) * sigma_m
J = gamma * [[0,0] [0,1.0+0im]]
He = [[deltaomega/2, 0.0] [0.0, 0.5*(-deltaomega - 1im*gamma) ]]

psi0 = zeros(ComplexF64, 2)
psi0[2] = 1 # Initial condition

ra_sys = System(H, # Hamiltonian
[sqrt(gamma)*sigma_m]) #Jump Operators
ra_params = SimulParameters(psi0,
    5.0, # Final time. Set very long so that all trajectories jump
    1, # seed
    1000, # Number of trajectories
    10_000, # Number of samples in the finegrid
    3, # Multiplier to use in the fine grid
    1e-3 # Tolerance for passing Dark state test
)
