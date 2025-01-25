import Distributions, HypothesisTests
using Test, BackAction, Statistics, LinearAlgebra
@testset verbose=true "Basic Operators" begin
       @test norm(BackAction.rd_sys.H - BackAction.rd_H) < BackAction.rd_EPS
       @test norm(BackAction.rd_sys.Ls[1]- BackAction.rd_L) < BackAction.rd_EPS
       @test norm(BackAction.rd_sys.Heff- BackAction.rd_He) < BackAction.rd_EPS
       @test norm(BackAction.rd_sys.J - BackAction.rd_J) < BackAction.rd_EPS
    end;

########### WTD ######################################
# Data generation
@time begin
rd_data = run_trajectories(BackAction.rd_sys, BackAction.rd_params)
end
rd_times = [rd_data[k][1].time for k in 1:BackAction.rd_params.ntraj if !isempty(rd_data[k])]
rd_d = Distributions.Exponential(1/BackAction.rd_gamma)
## Use a two sample Kolmogorov-Smirnov test, pvalue above 0.2 is accepted
rd_pvalue = HypothesisTests.pvalue(
HypothesisTests.ApproximateTwoSampleKSTest(rd_times, rand(rd_d, BackAction.rd_params.ntraj)))

@testset verbose=true "WTD (KS Test and fit)" begin
  rd_p0WTD = 0.2 # Minimal pvalue for accepting the null hypothesis
  fit_par = Distributions.fit(Distributions.Exponential, rd_times).θ
  @test rd_pvalue > rd_p0WTD
  @test abs(fit_par - 1/BackAction.rd_gamma) < 0.01
end

sys = BackAction.rd_sys
params = BackAction.rd_params
# Now from each trajectory, generate the states at the given times
sample_clicks = BackAction.run_trajectories(sys, params)
ntimes = 1000
t = collect(LinRange(0, params.tf, ntimes))
sample = Array{ComplexF64}(undef, sys.NLEVELS, ntimes, params.ntraj);
for n in 1:params.ntraj
    sample[:, :, n] = BackAction.states_att(t, sample_clicks[n], sys,  params.psi0)
end
# Check that all states are normalized
global local_flag = true
for n in 1:params.ntraj
    for k in 1:ntimes
        if (norm(sample[:, k, n]) - 1) > 0.01
            global local_flag = false
            break
        end
    end
end
# Check
accepted_error = 0.1 # Accept a relative accumulated error of up to 10%
# Obtain the observable on the sample
x_sample = zeros(ComplexF64, ntimes, params.ntraj)
for k in 1:params.ntraj
    for tn in 1:ntimes
        x_sample[tn, k] = dot(sample[:, tn, k], BackAction.sigma_z * sample[:, tn, k])
    end
end
x = real(dropdims( mean(x_sample, dims=2), dims=2));
x_theo = 2*exp.(-BackAction.rd_gamma.*t).-1
error = sum(abs.( (x - x_theo) ./ x_theo)) / (sum(abs.(x_theo)))

@test local_flag
@test error < accepted_error
