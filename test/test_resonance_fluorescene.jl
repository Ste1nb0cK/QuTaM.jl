
sys = QuTaM.rf_sys
params = QuTaM.rf_params
gamma = QuTaM.rf_gamma
delta = QuTaM.rf_delta
omega = QuTaM.rf_omega
r0 = [0.0; 0.0; -1.0] # Initial Condition
tspan = (0.0, params.tf)
t_given = collect(LinRange(0, params.tf, 1000));

# @testset "Resonance Fluorescene: Basic Operators" begin
#        @test norm(QuTaM.rf_sys.H - H) < QuTaM.rf_EPS
#        @test norm(QuTaM.rf_sys.Ls[1]- L) < QuTaM.rf_EPS
#        @test norm(QuTaM.rf_sys.Heff- Heff) < QuTaM.rf_EPS
#        @test norm(QuTaM.rf_sys.J - J) < QuTaM.rf_EPS

# end
@testset "Resonance Fluorescene: WTD" begin
################## Average Simulation ################3
# Generate a set of trajectories and states
sample_clicks = QuTaM.run_trajectories(sys, params)
ntimes = size(t_given)[1]
sample = zeros(ComplexF64, ntimes, sys.NLEVELS, params.ntraj) # states

for n in 1:params.ntraj
    states = QuTaM.evaluate_at_t(t_given, sample_clicks[n], sys,  params.psi0)
    for j in 1:sys.NLEVELS
        for tn in 1:ntimes
            sample[tn, j, n] = states[tn, j]
        end
    end
end

# Obtain the value of the observables
r_sample = zeros(Float64, ntimes, 3, params.ntraj)
sigma = [QuTaM.sigma_x, QuTaM.sigma_y, QuTaM.sigma_z]
for j in 1:params.ntraj
    for k in 1:3
        for tn in 1:ntimes
                r_sample[tn, k, j] = dot(sample[tn, :, j], sigma[k] * sample[tn, :, j])   # Drop the extra dimension
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
        # Replace with your custom formula
        gamma = d.gamma
        omega = d.omega
        return (16*gamma*omega^2)*exp(-0.5*gamma*tau) * sin(0.25*tau*sqrt(16*omega^2-gamma^2))^2/(-gamma^2+16*omega^2)
    end

function Distributions.cdf(d::WTD_rf, t::Real)
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

f = WTD_rf(gamma, omega) # Instance of the WTD
f_sample = rand(f, 500) # Sample from the WTD
# Run a KS test
pvalue = HypothesisTests.pvalue(HypothesisTests.ApproximateTwoSampleKSTest(tau_sample, f_sample))
rf_p0WTD = 0.2 # Minimal pvalue for accepting the null hypothesis
    @test pvalue > rf_p0WTD
end
