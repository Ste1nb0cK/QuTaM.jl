using DifferentialEquations

sys = QuTaM.rf_sys
params = QuTaM.rf_params
function rf_de!(dr, r, p, t)
    gamma = QuTaM.rf_gamma
    delta = QuTaM.rf_delta
    omega = QuTaM.rf_omega
    dr[1] = -0.5*gamma*r[1] - delta*r[2]
    dr[2] = delta*r[1] - 0.5*gamma*r[2] - omega*r[3]
    dr[3] = omega*r[2] - gamma*(r[3] + 1)
end

r0 = [0.0; 0.0; 1.0] # Initial Condition
tspan = (0.0, params.tf)
t_given = collect(LinRange(0, params.tf, 1000))
prob = ODEProblem(rf_de!, r0, tspan)
sol = solve(prob, reltol = 1e-6, saveat = t_given);

# Steady State
gamma = QuTaM.rf_gamma
delta = QuTaM.rf_delta
omega = QuTaM.rf_omega
r_steady = 1/(gamma^2 + 2*omega^2+4*delta^2) * [-4*delta*omega; 2*omega*gamma;-gamma^2-4*delta^2 ]
plot(sol, idxs=[(0,1), (0, 2), (0, 3)])
plot!(t_given, ones(Float64, size(t_given)[1])*r_steady[1], line=:dash)
plot!(t_given, ones(Float64, size(t_given)[1])*r_steady[2], line=:dash)
plot!(t_given, ones(Float64, size(t_given)[1])*r_steady[3], line=:dash)
