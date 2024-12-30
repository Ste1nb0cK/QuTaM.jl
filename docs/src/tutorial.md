# Tutorial
We will use a model of radiative damping of a two level atom as example.

## Basic Definitions
The basic way to use the library is to first define the system and 
the conditions of the simulation via the structs `System` and `SimulParameters` (see [System and SimulParameters](@ref basic_structs)).

### Defining the System
```@example workflow1
using QuTaM
# Define the system
deltaomega = 1
gamma = 1
H = 0.5*deltaomega * QuTaM.sigma_z # 2-level atom hamiltonian
L = sqrt(gamma) * QuTaM.sigma_m # Jump operator for Radiative Damping
sys = System(H, [L])
print(sys)
```

as you can see to create an instance of `System` one only needs the hamiltonian and the jump operators as it automatically computes all the derived quantities used in the simulation.

### Defining the Simulation

The key task the library has to do is sampling trajectories from the _Stochastic Master Equation_ (SME), from which the user will extract his 
statistics of interest. All the information necessary for generating this sample is contained in objects of the type  `SimulParameters`.

```@example workflow1
    
psi0 = zeros(ComplexF64, 2)
psi0[2] = 1

# Define the parameters
params = SimulParameters(psi0, # Initial condition
    3.0, # Final time. Set very long so that all trajectories jump
    1, # seed
    1000, # Number of trajectories
    50_000, # Number of samples in the finegrid
    10.5, # Multiplier to use in the fine grid
    1e-3 # Tolerance for passing Dark state test
)
print(params)
```
The two things that deserve commenting: 
1. The multiplier is used to generate the finegrid from which one samples the waiting times. Basically that finegrid is a discretization of `(0, params.multiplier * params.tf)` with step size `params.dt`.
2. In [Radaelli's paper](https://arxiv.org/abs/2303.15405) for the algorithm to fully work and have no problems in the statistics or infinite loops, one must have systems without darkstates. I don't think this is an absolutely necessart assumption, and it's probably possible to tweak it to work for dynamics with darksubspaces (as the one in this example), and so I added an improvised dark-state detection condition (to avoid infinite loops) which is passed if the probability of not jumping after a time `params.multiplier*params.tf` is inferior to the tolerance. See the source code of `run_single_trajectory` for details.

## Running the simulation

### ClickDetection and Trajectory
The data about a click is modeled by the struct `DetectionClick` [(more info)](@ref DetectionClick), this in practice works as a glorified tuple. Similarly, the library has an alias for `Vector{DetectionClick}` called `Trajectory`. In contrast with the paper, here the state after the jump is not used to fully identify the trajectory as the knowledge of the initial state is implied.

### Sampling trajectories
To obtain a sample of trajectories one uses `run_trajectories`, this returns a vector of trajectories as the sample :

```@example workflow1
data = run_trajectories(sys, params)
print(data[1:10])
```

We can access the clicks in a trajectory, **this is the main difference with things like QuTiP and QuantumOptics.jl**

```@example workflow1
traj = data[1] # Trajectory
click = traj[1] # first click in the trajectory
print(click.time,"\n") # time we waited to see the jump
print(click.label, "\n") # Channel in which the jump occured
```
## Obtaining the WTD
This is a simple example, and so we know in advance the WTD.

```@example workflow1
import Distributions, HypothesisTests # For the statistical test
using Plots

times = [data[k][1].time for k in 1:params.ntraj if !isempty(data[k])]
d = Distributions.Exponential(1/gamma) # Analytical result

# Plot
histogram(times, normalize=:pdf, label="Sample")
t = collect(LinRange(0, params.tf, 1000))
plot!(t, Distributions.pdf.(d, t), label="Analytical")
```

Looking good! to be completely sure we can do a [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test#Two-sample_Kolmogorov%E2%80%93Smirnov_test) and fit to an exponential.

```@example workflow1

HypothesisTests.ApproximateTwoSampleKSTest(times, rand(d, params.ntraj))
```
!!! note: In this type of tests the "sucess" corresponds with not rejecting the null hypothesis
```@example workflow1
fit_par = Distributions.fit(Distributions.Exponential, times).θ # Fit to an exponential
print("Estimated value:  $(fit_par) \n")
print("Real value: $(1/gamma) \n")
```

## Average 

Finally, we evaluate the expectation value of ``\sigma_z``. To achieve this we need to obtain the states in between jumps, this is done using `evaluate_at_t`. This function recieves a vector of given times and returns 
the states at the desired times.

```@example workflow1
using Statistics, LinearAlgebra
ntimes = 1000
sample = Array{ComplexF64}(undef, params.ntraj, ntimes, sys.NLEVELS ) # To store the states of all the trajectories
for n in 1:params.ntraj
        sample[n, :, :] = evaluate_at_t(t, data[n], sys,  params.psi0)
    end
# Obtain the observable on the sample
x_sample = zeros(ComplexF64, params.ntraj, ntimes)
for tn in 1:ntimes
    for k in 1:params.ntraj
        x_sample[k, tn] = dot(sample[k, tn, :], QuTaM.sigma_z * sample[k, tn, :])   # Drop the extra dimension
    end
end
x = real(dropdims( mean(x_sample, dims=1), dims=1));
x_theo = 2*exp.(-gamma.*t).-1

plot(t, x, label="Average")
plot!(t, x_theo, label="Analytical", line=:dash)
```
## Next steps
