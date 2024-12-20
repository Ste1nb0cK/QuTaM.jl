import Distributions, HypothesisTests
include("../src/Trajectories.jl")
using .QuTaM
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

sys = System(H, # Hamiltonian
[sqrt(gamma)*sigma_m]) #Jump Operators
simulparams = SimulParameters(psi0,
    5.0, # Final time. Set very long so that all trajectories jump
    1, # seed
    1000, # Number of trajectories
    10_000, # Number of samples in the finegrid
    3, # Multiplier to use in the fine grid
    1e-3 # Tolerance for passing Dark state test
)

@testset verbose=true "Basic Operators" begin
       @test norm(sys.H - H) < EPS
       @test norm(sys.Ls[1]- L) < EPS
       @test norm(sys.Heff- He) < EPS
       @test norm(sys.J - J ) < EPS
    end;
    # Data generation
data = run_trajectories(sys, simulparams)
# println(data[1][1].time)
times = [data[k][1].time for k in 1:simulparams.ntraj]
d = Distributions.Exponential(1/gamma)
## Use a two sample Kolmogorov-Smirnov test, pvalue above 0.2 is accepted
pvalue = HypothesisTests.pvalue(
   HypothesisTests.ApproximateTwoSampleKSTest(times, rand(d, simulparams.ntraj)))
@testset verbose=true "WTD Distribution" begin
   p0WTD = 0.2 # Minimal pvalue for accepting the null hypothesis
   @test pvalue > p0WTD
end
