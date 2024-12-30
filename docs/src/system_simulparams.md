# [System and SimulParameters](@id basic_structs)

The two basic parts of any workflow are the structs `System` and `SimulParameters`. 
The former defines the dynamics and the latter the specifics of the simulation e.g.
grid precision and number of trajectories.

```@docs; canonical=false
System
```

```@docs; canonical=false
SimulParameters
```

# Example
Below we give an example using the predefined Pauli matrices available in the package.

```@example
using QuTaM

# Define the system
deltaomega = 1
gamma = 1
H = 0.5*deltaomega * QuTaM.sigma_z # 2-level atom hamiltonian
L = sqrt(gamma) * QuTaM.sigma_m # Jump operator for Radiative Damping
sys = System(H, [L])

psi0 = zeros(ComplexF64, 2)
psi0[2] = 1 # Initial condition

# Define the parameters

params = SimulParameters(psi0,
    3.0, # Final time. Set very long so that all trajectories jump
    1, # seed
    1000, # Number of trajectories
    50_000, # Number of samples in the finegrid
    10.5, # Multiplier to use in the fine grid
    1e-3 # Tolerance for passing Dark state test
)

print(sys, "\n\n")
print(params, "\n")
```

# [The DetectionClick Struct](@id  )

To implement the clicks, the struct `DetectionClick` is used.

```@docs; canonical=false
DetectionClick
```
```@docs; canonical=false
Trajectory
```


