import HypothesisTests, Distributions
using QuTaM, LinearAlgebra, Statistics, Random, QuadGK, Test, Plots,
    OrdinaryDiffEq
########## Setup
sys = QuTaM.rf_sys
params = QuTaM.rf_params
r0 = [0.0; 0.0; -1.0] # Initial Condition
tspan = (0.0, params.tf)
t_given = collect(LinRange(0, params.tf, 1000));

################## Average Simulation ################3
# Generate a set of trajectories and states
println("Sampling clicks\n")
@time begin
    sample_clicks = QuTaM.run_trajectories(sys, params)
end
ntimes = size(t_given)[1]
sample = zeros(ComplexF64,  sys.NLEVELS, ntimes, params.ntraj) # states

println("Obtaining states between jumps\n")
@time begin
    for n in 1:params.ntraj
        states = QuTaM.evaluate_at_t(t_given, sample_clicks[n], sys,  params.psi0)
                sample[:, :, n] = states[:, :]
    end
end
# Obtain the value of the observables
r_sample = zeros(Float64, ntimes, 3, params.ntraj)
sigma = [QuTaM.sigma_x, QuTaM.sigma_y, QuTaM.sigma_z]
println("Calculating expectation of observables\n")

@time begin
    for j in 1:params.ntraj
        for k in 1:3
            for tn in 1:ntimes
                r_sample[tn, k, j] = dot(sample[:, tn, j], sigma[k] * sample[:, tn, j])   # Drop the extra dimension
            end
        end
    end
end
r_avg = dropdims(mean(r_sample, dims=3), dims=3) # Ensemble average

# Exploiting that this is a renewal process, obtain a sample of waiting times
tau_sample = Vector{Float64}()
for traj in sample_clicks
    if !isempty(traj)
        for click in traj
            push!(tau_sample, click.time)
        end
    else
        continue
    end
end

# Define the analytical result for the WTD
struct WTD_rf <: Distributions.ContinuousUnivariateDistribution
    omega::Float64
    gamma::Float64

end

function Distributions.support(d::WTD_rf)
    return Distributions.Interval(0, Inf)
end

function Distributions.pdf(d::WTD_rf, tau::Real)
        gamma = d.gamma
        omega = d.omega
        return (16*gamma*omega^2)*exp(-0.5*gamma*tau) * sin(0.25*tau*sqrt(16*omega^2-gamma^2))^2/(-gamma^2+16*omega^2)
    end

function Distributions.cdf(d::WTD_rf, t::Real)
    gamma = d.gamma
    omega = d.omega
    pdf(tau) = (16*gamma*omega^2)*exp(-0.5*gamma*tau) * sin(0.25*tau*sqrt(16*omega^2-gamma^2))^2/(-gamma^2+16*omega^2)
    return quadgk(pdf, 0, t, rtol=1e-8)[1]
end

function Base.rand(rng::AbstractRNG, d::WTD_rf)
    # Use inversion sampling
    alpha = rand()
    t = 0
    dt = 0.001 # This is a magic number, the point is that this matches the dt in rf_params
    while Distributions.cdf(d, t) < alpha
        t = dt + t
    end
    return t  # Return a sample
end

f = WTD_rf(QuTaM.rf_gamma, QuTaM.rf_omega) # Instance of the WTD
println("Sampling from the analytical WTD\n")
@time begin
    f_sample = rand(f, 500) # Sample from the WTD
end
# Run a KS test
println("Runnign KS test\n")
@time begin
    pvalue = HypothesisTests.pvalue(HypothesisTests.ApproximateTwoSampleKSTest(tau_sample, f_sample))
end
rf_p0WTD = 0.2 # Minimal pvalue for accepting the null hypothesis

@testset verbose=true "WTD" begin
    @test pvalue > rf_p0WTD
end

################# Visual test of the WTD #############
hist = histogram(tau_sample, normalize=:pdf, label="Sample")
plot!(t_given, Distributions.pdf.(f, t_given), label="Analytical")
################ Visual test of the observables
# System of equations for RF (Source: Wiseman section 3.3.1)
function rf_de!(dr, r, p, t)
    gamma = QuTaM.rf_gamma
    delta = QuTaM.rf_delta
    omega = QuTaM.rf_omega
    dr[1] = -0.5*gamma*r[1] - 2*delta*r[2]
    dr[2] = 2*delta*r[1] - 0.5*gamma*r[2] - 2*omega*r[3]
    dr[3] = 2*omega*r[2] - gamma*(r[3] + 1)
end

# Solution to the unconditional equation
prob = ODEProblem(rf_de!, r0, tspan)
sol = solve(prob, reltol = 1e-6, saveat = t_given);

using LaTeXStrings

p1 = plot(t_given, r_avg[:, 1], label="Unraveling", title=L"\sigma_x", seriescolor=:blue)
plot!(t_given, r_sample[:,1, 1], label="Sample Trajectory", seriescolor=:green)
plot!(sol, idxs =(0, 1), label="Lindblad", seriescolor=:red)

p2 = plot(t_given, r_avg[:, 2], title=L"\sigma_y", seriescolor=:blue)
plot!(t_given, r_sample[:,2, 1], seriescolor=:green)
plot!(sol, idxs =(0, 2),  legend=false,  seriescolor=:red)

p3 = plot(t_given, r_avg[:, 3], title=L"\sigma_z", seriescolor=:blue)
plot!(t_given, r_sample[:,3, 1],  seriescolor=:green)
plot!(sol, idxs =(0, 3), legend=false, seriescolor=:red)

l = @layout [
    a{0.5w} grid(3,1)
]

plot(hist, p1, p2,p3,  layout=l)
