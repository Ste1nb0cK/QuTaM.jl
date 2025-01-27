# Exported Functions
```@docs
monitoringoperator(t_given::Vector{Float64},
    sys::System, Heff_par::Function, Ls_par, traj::Trajectory, psi0::Vector{ComplexF64}, theta::Vector{Float64},
                            dtheta::Vector{Float64})
```

```@docs
BackAction.getheff_parametrized(H_par::Function, Ls_par)
```

## Internal Functions

```@docs

BackAction.expheff_derivative(Heff_par::Function, tau::Float64, theta::Vector{Float64}, dtheta::Vector{Float64})
```

```@docs
BackAction.jumpoperators_derivatives(Ls_par, theta, dtheta)
```

```@docs
BackAction.writederivative!(dpsi::SubArray{ComplexF64, 1}, L::Matrix{ComplexF64}, dL::SubArray{ComplexF64, 2}, V::Matrix{ComplexF64}, dV::Matrix{ComplexF64},
psi0::Vector{ComplexF64})
```

```@docs
BackAction.writederivative!(dpsi::SubArray{ComplexF64, 1},L::Matrix{ComplexF64},dL::SubArray{ComplexF64, 2},V::Matrix{ComplexF64}, dV::Matrix{ComplexF64},psi0::SubArray{ComplexF64, 1},dpsi0::SubArray{ComplexF64, 1})
```

```@docs
BackAction.derivatives_atjumps(sys::System, Heff_par::Function, Ls_par, traj::Trajectory, psi0::Vector{ComplexF64}, theta::Vector{Float64}, dtheta::Vector{Float64})
```

```@docs
BackAction.writexi!(xi::SubArray{ComplexF64, 2}, dV::Matrix{ComplexF64},psi::SubArray{ComplexF64, 1}, psi0::Vector{ComplexF64})
```

```@docs
BackAction.writexi!(xi::SubArray{ComplexF64, 2}, V::Matrix{ComplexF64}, dV::Matrix{ComplexF64}, psijump::SubArray{ComplexF64, 1}, dpsijump::SubArray{ComplexF64, 1},psi::SubArray{ComplexF64, 1})
```
