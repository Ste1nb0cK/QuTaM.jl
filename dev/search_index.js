var documenterSearchIndex = {"docs":
[{"location":"system_simulparams.html#basic_structs","page":"System and SimulParameters","title":"System and SimulParameters","text":"","category":"section"},{"location":"system_simulparams.html","page":"System and SimulParameters","title":"System and SimulParameters","text":"The two basic parts of any workflow are the structs System and SimulParameters.  The former defines the dynamics and the latter the specifics of the simulation e.g. grid precision and number of trajectories.","category":"page"},{"location":"system_simulparams.html","page":"System and SimulParameters","title":"System and SimulParameters","text":"System","category":"page"},{"location":"system_simulparams.html#QuTaM.System-system_simulparams","page":"System and SimulParameters","title":"QuTaM.System","text":"System(\nNLEVELS::Int64, NCHANNELS::Int64, H::Matrix{ComplexF64}\nLs::Vector{Matrix{ComplexF64}}, J::Matrix{ComplexF64},\nHeff::Matrix{ComplexF64})\n\nA mutable struct that characterizes the dynamics via specification of  the jump and hamiltonian operators.\n\nFields\n\nNLEVELS::Int64: Number of levels of the system\nNCHANNELS::Int64: Number of jump channels\nH::Matrix{ComplexF64}: Hamiltonian\nLs::Vector{Matrix{ComplexF64}}: List of jump operators\nJ::Matrix{ComplexF64}: Sum of all the L_k^*L_k\nHeff::Matrix{ComplexF64}: Effective Hamiltonian\n\nConstructor\n\nTo create an instance it's enough to provide the hamiltonian and the jump operators in a vector. System(H::Matrix{ComplexF64}, Ls::Vector{Matrix{ComplexF64}})\n\n\n\n\n\n","category":"type"},{"location":"system_simulparams.html","page":"System and SimulParameters","title":"System and SimulParameters","text":"SimulParameters","category":"page"},{"location":"system_simulparams.html#QuTaM.SimulParameters-system_simulparams","page":"System and SimulParameters","title":"QuTaM.SimulParameters","text":"SimulParameters(\n    psi0::Vector{ComplexF64}, nsamples::Int64, seed::Int64,\n            ntraj::Int64, multiplier::Float64, tf::Float64,\n            dt::Float64, eps::Float64)\n\nA mutable struct containing all the necessary information for running the the simulation.\n\nFields\n\npsi0::Vector{ComplexF64}: Initial state vector\nnsamples::Int64: Number of samples in the finegrid\nseed::Int64: seed\nntraj::Int64: Number of trajectories\nmultiplier::Float64: Multiplier to use in the fine grid\ntf::Float64: Final time\ndt::Float64: time step for the finegrid\neps::Float64: Tolerance for passing WTD normalziation\n\nConstructor\n\nTo create an instance it's enough to provide initial state, final time, seed and number of trajectories. Unless given nsamples, multiplier and eps use default values. SimulParameters(psi0::Vector{ComplexF64}, tf::Float64,         s::Int64, ntraj::Int64, nsamples::Int64=10000, m::Float64=10.0,                              eps::Float64=1e-3)\n\nAbout the multiplier\n\nFor the Gillipsie algorithm to work it's key to have a grid that's capable of resolving the statistical details of the WTD, this grid is taken in the interval (0, tf*multiplier).\n\n\n\n\n\n","category":"type"},{"location":"system_simulparams.html#Example","page":"System and SimulParameters","title":"Example","text":"","category":"section"},{"location":"system_simulparams.html","page":"System and SimulParameters","title":"System and SimulParameters","text":"Below we give an example using the predefined Pauli matrices available in the package.","category":"page"},{"location":"system_simulparams.html","page":"System and SimulParameters","title":"System and SimulParameters","text":"using QuTaM\n\n# Define the system\ndeltaomega = 1\ngamma = 1\nH = 0.5*deltaomega * QuTaM.sigma_z # 2-level atom hamiltonian\nL = sqrt(gamma) * QuTaM.sigma_m # Jump operator for Radiative Damping\nsys = System(H, [L])\n\npsi0 = zeros(ComplexF64, 2)\npsi0[2] = 1 # Initial condition\n\n# Define the parameters\n\nparams = SimulParameters(psi0,\n    3.0, # Final time. Set very long so that all trajectories jump\n    1, # seed\n    1000, # Number of trajectories\n    50_000, # Number of samples in the finegrid\n    10.5, # Multiplier to use in the fine grid\n    1e-3 # Tolerance for passing Dark state test\n)\n\nprint(sys, \"\\n\\n\")\nprint(params, \"\\n\")","category":"page"},{"location":"system_simulparams.html#","page":"System and SimulParameters","title":"The DetectionClick Struct","text":"","category":"section"},{"location":"system_simulparams.html","page":"System and SimulParameters","title":"System and SimulParameters","text":"To implement the clicks, the struct DetectionClick is used.","category":"page"},{"location":"system_simulparams.html","page":"System and SimulParameters","title":"System and SimulParameters","text":"DetectionClick","category":"page"},{"location":"system_simulparams.html#QuTaM.DetectionClick-system_simulparams","page":"System and SimulParameters","title":"QuTaM.DetectionClick","text":"DetectionClick(time::Float64, label::Int64)\n\nInmutable struct that represents the clicks by the time waited to see the click and the label of the channel in which it occured.\n\nFields\n\ntime::Float64: Waiting time\nlabel::Int64: Label of the channel of the click\n\n\n\n\n\n\n\n","category":"type"},{"location":"system_simulparams.html","page":"System and SimulParameters","title":"System and SimulParameters","text":"Trajectory","category":"page"},{"location":"system_simulparams.html#QuTaM.Trajectory-system_simulparams","page":"System and SimulParameters","title":"QuTaM.Trajectory","text":"Alias for Vector{DetectionClick}\n\n\n\n\n\n","category":"type"},{"location":"reference.html#reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference.html","page":"Reference","title":"Reference","text":"Here you can find the documentation for the all the functions and structs in the library.","category":"page"},{"location":"reference.html","page":"Reference","title":"Reference","text":"Modules = [QuTaM]","category":"page"},{"location":"reference.html#QuTaM.DetectionClick","page":"Reference","title":"QuTaM.DetectionClick","text":"DetectionClick(time::Float64, label::Int64)\n\nInmutable struct that represents the clicks by the time waited to see the click and the label of the channel in which it occured.\n\nFields\n\ntime::Float64: Waiting time\nlabel::Int64: Label of the channel of the click\n\n\n\n\n\n\n\n","category":"type"},{"location":"reference.html#QuTaM.SimulParameters","page":"Reference","title":"QuTaM.SimulParameters","text":"Inner constructor of SimulParameters SimulParameters(psi0::Vector{ComplexF64}, tf::Float64,         s::Int64, ntraj::Int64, nsamples::Int64=10000, m::Float64=10.0,                              eps::Float64=1e-3)\n\n\n\n\n\n","category":"type"},{"location":"reference.html#QuTaM.SimulParameters-2","page":"Reference","title":"QuTaM.SimulParameters","text":"SimulParameters(\n    psi0::Vector{ComplexF64}, nsamples::Int64, seed::Int64,\n            ntraj::Int64, multiplier::Float64, tf::Float64,\n            dt::Float64, eps::Float64)\n\nA mutable struct containing all the necessary information for running the the simulation.\n\nFields\n\npsi0::Vector{ComplexF64}: Initial state vector\nnsamples::Int64: Number of samples in the finegrid\nseed::Int64: seed\nntraj::Int64: Number of trajectories\nmultiplier::Float64: Multiplier to use in the fine grid\ntf::Float64: Final time\ndt::Float64: time step for the finegrid\neps::Float64: Tolerance for passing WTD normalziation\n\nConstructor\n\nTo create an instance it's enough to provide initial state, final time, seed and number of trajectories. Unless given nsamples, multiplier and eps use default values. SimulParameters(psi0::Vector{ComplexF64}, tf::Float64,         s::Int64, ntraj::Int64, nsamples::Int64=10000, m::Float64=10.0,                              eps::Float64=1e-3)\n\nAbout the multiplier\n\nFor the Gillipsie algorithm to work it's key to have a grid that's capable of resolving the statistical details of the WTD, this grid is taken in the interval (0, tf*multiplier).\n\n\n\n\n\n","category":"type"},{"location":"reference.html#QuTaM.System","page":"Reference","title":"QuTaM.System","text":"System(\nNLEVELS::Int64, NCHANNELS::Int64, H::Matrix{ComplexF64}\nLs::Vector{Matrix{ComplexF64}}, J::Matrix{ComplexF64},\nHeff::Matrix{ComplexF64})\n\nA mutable struct that characterizes the dynamics via specification of  the jump and hamiltonian operators.\n\nFields\n\nNLEVELS::Int64: Number of levels of the system\nNCHANNELS::Int64: Number of jump channels\nH::Matrix{ComplexF64}: Hamiltonian\nLs::Vector{Matrix{ComplexF64}}: List of jump operators\nJ::Matrix{ComplexF64}: Sum of all the L_k^*L_k\nHeff::Matrix{ComplexF64}: Effective Hamiltonian\n\nConstructor\n\nTo create an instance it's enough to provide the hamiltonian and the jump operators in a vector. System(H::Matrix{ComplexF64}, Ls::Vector{Matrix{ComplexF64}})\n\n\n\n\n\n","category":"type"},{"location":"reference.html#QuTaM.System-Tuple{Matrix{ComplexF64}, Vector{Matrix{ComplexF64}}}","page":"Reference","title":"QuTaM.System","text":"     Inner Constructor of `System` struct.\n     # Arguments:\n     `H::Matrix{ComplexF64}`\n     `Ls::Vector{Matrix{ComplexF64}}`\n\n\n\n\n\n","category":"method"},{"location":"reference.html#QuTaM.Trajectory","page":"Reference","title":"QuTaM.Trajectory","text":"Alias for Vector{DetectionClick}\n\n\n\n\n\n","category":"type"},{"location":"reference.html#QuTaM.precompute!-Tuple{System, Int64, Vector{Float64}, Vector{Matrix{ComplexF64}}}","page":"Reference","title":"QuTaM.precompute!","text":"precompute!(sys::System, nsamples::Int64,\n     ts::Vector{Float64}, Qs::Vector{Matrix{ComplexF64}})\n\nDoes the precomputation routine for the Gillipsie algorithm. The result is stored in the Qs.\n\nArguments\n\nsys::System: system nsamples::Int64: number of samples ts::Vector{Float64}: fine grid vector. IT IS ASSUMED HOMOGENOUS Qs::Vector{Matrix{ComplexF64}}: vector of matrices to store the precomputation.\n\n\n\n\n\n","category":"method"},{"location":"reference.html#QuTaM.run_single_trajectory-Tuple{System, SimulParameters, Vector{Float64}, Vector{Float64}, Vector{ComplexF64}, Vector{Float64}, Vector{Matrix{ComplexF64}}}","page":"Reference","title":"QuTaM.run_single_trajectory","text":"run_single_trajectories(sys::System, params::SimulParameters,\n W::Vector{Float64}, P::Vector{Float64}, psi::Vector{ComplexF64},\n ts::Vector{Float64}, Qs::Vector{Matrix{ComplexF64}};\n  seed::Int64 = 1) -> Trajectory\n\nSample a single trajectory from the system and parameters using the Gillipsie algorithm . This is inteded to be used by run_trajectories\n\nArguments:\n\nsys::System: System of interest\nparams::SimulParameters: simulation parameters\nW::Vector{Float64}: to store the weights over the fine grid. Its lenght must be 1 more than that of ts\nP::Vector{Float64}: to store the weights over the channels\npsi::Vector{ComplexF64}: to store the current state vector\nts::Vector{Float64}: the finegrid of waiting times\nQs::Vector{Matrix{ComplexF64}}: to store the precomputed values\nseed::Int64: seed for generating the trajectory\n\nReturns:\n\nThe sample trajectory\n\nWarning: the seed does not coincide with that of params by default.\n\nWaning: final jump in the trajectory happens after final time\n\nThe trajectory ends when a jump that happens after params.tf is obtained, yet that jump is stored in the trajectory. In other words, the last jump of the trajectory always happen after the set final time.\n\n\n\n\n\n","category":"method"},{"location":"reference.html#QuTaM.run_trajectories-Tuple{System, SimulParameters}","page":"Reference","title":"QuTaM.run_trajectories","text":"run_trajectories(sys::System, params::SimulParameters) -> Vector{Trajectory}\n\nSample multiple trajectories for a given system and parameters.\n\nArguments\n\nsys::System: The quantum system to simulate, containing information about its structure, energy levels, and dynamics.\nparams::SimulParameters: A structure containing simulation parameters such as:\nprogbar::Bool: show progress bar or not. true by default.\n\nReturns\n\nVector{Trajectory}: A vector containing the results of the simulated trajectories. Each element corresponds to a single trajectory and encapsulates relevant system state information over time.\n\n\n\n\n\n","category":"method"},{"location":"reference.html#QuTaM.sample_single_trajectory-Tuple{System, SimulParameters, Int64}","page":"Reference","title":"QuTaM.sample_single_trajectory","text":"sample_single_trajectory(sys::System, params::SimulParameters, seed::Int64) -> Vector{Trajectory}\n\nSample a single trajectory. This is inteded to be used in case a single sample is needed without redefining params.ntraj.\n\nArguments\n\nsys::System\nparams::SimulParameters\nseed::Int64\n\nReturnrs\n\nA 1-element trajectory vector.\n\n\n\n\n\n","category":"method"},{"location":"reference.html#QuTaM.states_at_jumps-Tuple{Vector{DetectionClick}, System, Vector{ComplexF64}}","page":"Reference","title":"QuTaM.states_at_jumps","text":"states_at_jumps(traj::Trajectory, sys::System,\n                  psi0::Vector{ComplexF64}) - > Array{ComplexF64}\n\nFrom a given trajectory, recover the states at each jump.\n\nArguments\n\ntraj::Trajectory: Trajectory\nsys::System: System of interest\npsi0::Vector{ComplexF64}: Initial state\n\nReturns\n\nArray{ComplexF64} with the states\n\nThe dimensions of the returned array s are (size(traj), sys.NLEVELS), so to recover the state vector at the n-th jump one would do s[n, :].\n\n\n\n\n\n","category":"method"},{"location":"index.html","page":"Home","title":"Home","text":"CurrentModule = QuTaM","category":"page"},{"location":"index.html#QuTaM","page":"Home","title":"QuTaM","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Here I would write a good description of the library, if I had one.","category":"page"},{"location":"tutorial.html#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"We will use a model of radiative damping of a two level atom as example.","category":"page"},{"location":"tutorial.html#Basic-Definitions","page":"Tutorial","title":"Basic Definitions","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"The basic way to use the library is to first define the system and  the conditions of the simulation via the structs System and SimulParameters (see System and SimulParameters).","category":"page"},{"location":"tutorial.html#Defining-the-System","page":"Tutorial","title":"Defining the System","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"using QuTaM\n# Define the system\ndeltaomega = 1\ngamma = 1\nH = 0.5*deltaomega * QuTaM.sigma_z # 2-level atom hamiltonian\nL = sqrt(gamma) * QuTaM.sigma_m # Jump operator for Radiative Damping\nsys = System(H, [L])\nprint(sys)","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"as you can see to create an instance of System one only needs the hamiltonian and the jump operators as it automatically computes all the derived quantities used in the simulation.","category":"page"},{"location":"tutorial.html#Defining-the-Simulation","page":"Tutorial","title":"Defining the Simulation","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"The key task the library has to do is sampling trajectories from the Stochastic Master Equation (SME), from which the user will extract his  statistics of interest. All the information necessary for generating this sample is contained in objects of the type  SimulParameters.","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"    \npsi0 = zeros(ComplexF64, 2)\npsi0[2] = 1\n\n# Define the parameters\nparams = SimulParameters(psi0, # Initial condition\n    3.0, # Final time. Set very long so that all trajectories jump\n    1, # seed\n    1000, # Number of trajectories\n    50_000, # Number of samples in the finegrid\n    10.5, # Multiplier to use in the fine grid\n    1e-3 # Tolerance for passing Dark state test\n)\nprint(params)","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"The two things that deserve commenting: ","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"The multiplier is used to generate the finegrid from which one samples the waiting times. Basically that finegrid is a discretization of (0, params.multiplier * params.tf) with step size params.dt.\nIn Radaelli's paper for the algorithm to fully work and have no problems in the statistics or infinite loops, one must have systems without darkstates. I don't think this is an absolutely necessart assumption, and it's probably possible to tweak it to work for dynamics with darksubspaces (as the one in this example), and so I added an improvised dark-state detection condition (to avoid infinite loops) which is passed if the probability of not jumping after a time params.multiplier*params.tf is inferior to the tolerance. See the source code of run_single_trajectory for details.","category":"page"},{"location":"tutorial.html#Running-the-simulation","page":"Tutorial","title":"Running the simulation","text":"","category":"section"},{"location":"tutorial.html#ClickDetection-and-Trajectory","page":"Tutorial","title":"ClickDetection and Trajectory","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"The data about a click is modeled by the struct DetectionClick (more info), this in practice works as a glorified tuple. Similarly, the library has an alias for Vector{DetectionClick} called Trajectory. In contrast with the paper, here the state after the jump is not used to fully identify the trajectory as the knowledge of the initial state is implied.","category":"page"},{"location":"tutorial.html#Sampling-trajectories","page":"Tutorial","title":"Sampling trajectories","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"To obtain a sample of trajectories one uses run_trajectories, this returns a vector of trajectories as the sample :","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"data = run_trajectories(sys, params)\nprint(data[1:10])","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"We can access the clicks in a trajectory, this is the main difference with things like QuTiP and QuantumOptics.jl","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"traj = data[1] # Trajectory\nclick = traj[1] # first click in the trajectory\nprint(click.time,\"\\n\") # time we waited to see the jump\nprint(click.label, \"\\n\") # Channel in which the jump occured","category":"page"},{"location":"tutorial.html#Obtaining-the-WTD","page":"Tutorial","title":"Obtaining the WTD","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"This is a simple example, and so we know in advance the WTD.","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"import Distributions, HypothesisTests # For the statistical test\nusing Plots\n\ntimes = [data[k][1].time for k in 1:params.ntraj if !isempty(data[k])]\nd = Distributions.Exponential(1/gamma) # Analytical result\n\n# Plot\nhistogram(times, normalize=:pdf, label=\"Sample\")\nt = collect(LinRange(0, params.tf, 1000))\nplot!(t, Distributions.pdf.(d, t), label=\"Analytical\")","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"Looking good! to be completely sure we can do a Kolmogorov-Smirnov test and fit to an exponential.","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"\nHypothesisTests.ApproximateTwoSampleKSTest(times, rand(d, params.ntraj))","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"!!! note: In this type of tests the \"sucess\" corresponds with not rejecting the null hypothesis","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"fit_par = Distributions.fit(Distributions.Exponential, times).θ # Fit to an exponential\nprint(\"Estimated value:  $(fit_par) \\n\")\nprint(\"Real value: $(1/gamma) \\n\")","category":"page"},{"location":"tutorial.html#Average","page":"Tutorial","title":"Average","text":"","category":"section"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"Finally, we evaluate the expectation value of sigma_z. To achieve this we need to obtain the states in between jumps, this is done using evaluate_at_t. This function recieves a vector of given times and returns  the states at the desired times.","category":"page"},{"location":"tutorial.html","page":"Tutorial","title":"Tutorial","text":"using Statistics, LinearAlgebra\nntimes = 1000\nsample = Array{ComplexF64}(undef, params.ntraj, ntimes, sys.NLEVELS ) # To store the states of all the trajectories\nfor n in 1:params.ntraj\n        sample[n, :, :] = evaluate_at_t(t, data[n], sys,  params.psi0)\n    end\n# Obtain the observable on the sample\nx_sample = zeros(ComplexF64, params.ntraj, ntimes)\nfor tn in 1:ntimes\n    for k in 1:params.ntraj\n        x_sample[k, tn] = dot(sample[k, tn, :], QuTaM.sigma_z * sample[k, tn, :])   # Drop the extra dimension\n    end\nend\nx = real(dropdims( mean(x_sample, dims=1), dims=1));\nx_theo = 2*exp.(-gamma.*t).-1\n\nplot(t, x, label=\"Average\")\nplot!(t, x_theo, label=\"Analytical\", line=:dash)","category":"page"},{"location":"tutorial.html#Next-steps","page":"Tutorial","title":"Next steps","text":"","category":"section"}]
}