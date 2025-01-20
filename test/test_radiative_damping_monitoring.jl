using QuTaM, Plots, Test, LinearAlgebra, Statistics
sys = QuTaM.rd_sys
params = SimulParameters(QuTaM.rd_psi0,
    3.0, # Final time. Set very long so that all trajectories jump
    1, # seed
    1000, # Number of trajectories
    50_000, # Number of samples in the finegrid
    10.5, # Multiplier to use in the fine grid
    1e-3 # Tolerance for passing Dark state test
)

trajectories = run_trajectories(sys, params);
# Parametrization stuff
H_parametrized = (delta::Float64, gamma::Float64) -> (0.5*delta*QuTaM.sigma_z)::Matrix{ComplexF64}
L_parametrized = (delta::Float64, gamma::Float64) -> (sqrt(gamma)*QuTaM.sigma_m)::Matrix{ComplexF64}
Heff_parametrized = GetHeffParametrized(H_parametrized, [L_parametrized])

ntimes = 100
t_given = collect(LinRange(0, params.tf, ntimes));

# Obtain Monitoring Operator
xi_sample = Array{ComplexF64}(undef, sys.NLEVELS, sys.NLEVELS, ntimes, params.ntraj)
for n in 1:params.ntraj
    xi_sample[:, :, :, n] = MonitoringOperator(t_given, sys, Heff_parametrized, [L_parametrized], trajectories[n], params.psi0,
                           [QuTaM.rd_deltaomega, QuTaM.rd_gamma], [0.0, QuTaM.rd_gamma/100])
end

# Calculate the sample Fisher Information
fi_sample = Array{Float64}(undef, ntimes, params.ntraj)
for n in 1:params.ntraj
    for k in 1:ntimes
       fi_sample[k, n] = real(tr(xi_sample[:, :, k, n]))^2
    end
end
# Fisher Information Average
fi = dropdims(mean(fi_sample, dims=2), dims=2);

# Check global error against analytical result
f_analytical(t) = (1-exp(-QuTaM.rd_gamma*t))/(QuTaM.rd_gamma^2)
fi_theo = f_analytical.(t_given)

fi_min, fi_max = extrema(fi_theo)
# Use as error measure the normalized MSRE
accepted_error = 0.1
error_global = sqrt( mean((fi_theo - fi).^2 )) /(fi_max - fi_min)

@test error_global < accepted_error

# Plot against analytical result
plot(t_given, f_analytical.(t_given), label="Analytical", color="black", linewidth=3.3)
scatter!(t_given, fi, label="simulation", color="blue")
