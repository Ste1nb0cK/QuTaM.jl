import Distributions, HypothesisTests
include("../src/Trajectories.jl")
using .QuTaM

@testset verbose=true "Basic Operators" begin
       @test norm(QuTaM.rd_sys.H - QuTaM.rd_H) < QuTaM.rd_EPS
       @test norm(QuTaM.rd_sys.Ls[1]- QuTaM.rd_L) < QuTaM.rd_EPS
       @test norm(QuTaM.rd_sys.Heff- QuTaM.rd_He) < QuTaM.rd_EPS
       @test norm(QuTaM.rd_sys.J - QuTaM.rd_J) < QuTaM.rd_EPS
    end;
    # Data generation
rd_data = run_trajectories(QuTaM.rd_sys, QuTaM.rd_params)
# println(data[1][1].time)
rd_times = [rd_data[k][1].time for k in 1:QuTaM.rd_params.ntraj]
rd_d = Distributions.Exponential(1/QuTaM.rd_gamma)
## Use a two sample Kolmogorov-Smirnov test, pvalue above 0.2 is accepted
rd_pvalue = HypothesisTests.pvalue(
   HypothesisTests.ApproximateTwoSampleKSTest(rd_times, rand(rd_d, QuTaM.rd_params.ntraj)))
@testset verbose=true "WTD Distribution" begin
   rd_p0WTD = 0.2 # Minimal pvalue for accepting the null hypothesis
   @test rd_pvalue > rd_p0WTD
   println("WTD Distribution fitted with pvalue=$(rd_pvalue)\n")
end
