function getheff_parametrized(H_par::Function, Ls_par)::Function
    # here theta is expected to be a vector
    return (theta...) -> begin # Get an arbitrary number of arguments
        # The ... "splatts" the vector so it is passed as a tuple to the function
        LLs_par = [ adjoint(L_par(theta...))*L_par(theta...) for L_par in Ls_par]
        return H_par(theta...) - 0.5im*sum(LLs_par)
    end
end

function expheff_derivative(Heff_par::Function, tau::Float64, theta::Vector{Float64}, dtheta::Vector{Float64})
    f1 =  exp(-1im*tau*Heff_par((theta + 2*dtheta)...))
    f2 =  exp(-1im*tau*Heff_par((theta + 1*dtheta)...))
    f3 =  exp(-1im*tau*Heff_par((theta - 1*dtheta)...))
    f4 =  exp(-1im*tau*Heff_par((theta - 2*dtheta)...))
    return (-f1 + 8*f2 - 8*f3 + f4 )/(12*norm(dtheta))
end

function jumpoperators_derivatives(Ls_par, theta, dtheta)
    nchannels = size(Ls_par)[1]
    nlevels = size(Ls_par[1](theta...))[1]
    dLs = zeros(ComplexF64, nlevels, nlevels, nchannels)
    for k in 1:nchannels
        f1 = Ls_par[k]((theta + 2*dtheta)...)
        f2 = Ls_par[k]((theta+dtheta)...)
        f3 = Ls_par[k]((theta-dtheta)...)
        f4 = Ls_par[k]((theta-2*dtheta)...)
        dLs[:, :, k] = (-f1 + 8*f2 - 8*f3 + f4 )/(12*norm(dtheta))
    end
    return dLs
end

function writederivative!(dpsi::SubArray{ComplexF64, 1},
                          L::Matrix{ComplexF64},
                          dL::SubArray{ComplexF64, 2},
                          V::Matrix{ComplexF64}, dV::Matrix{ComplexF64},
                          psi0::Vector{ComplexF64})
    dpsi .= (dL*V + L*dV)*psi0
end

function writederivative!(dpsi::SubArray{ComplexF64, 1},
                          L::Matrix{ComplexF64},
                          dL::SubArray{ComplexF64, 2},
                          V::Matrix{ComplexF64}, dV::Matrix{ComplexF64},
                          psi0::SubArray{ComplexF64, 1},#Union{Vector{ComplexF64}, SubArray{ComplexF64}},
                          dpsi0::SubArray{ComplexF64, 1})
    dpsi .= (dL*V + L*dV)*psi0 + L*V*dpsi0
end

function derivatives_atjumps(sys::System, Heff_par::Function, Ls_par, traj::Trajectory, psi0::Vector{ComplexF64}, theta::Vector{Float64},
                            dtheta::Vector{Float64})
    # 0. Special Case: if the trajectory is empty, return an empty array
    if isempty(traj)
        return Array{ComplexF64}(undef, 0, 0)
    end
    # 1. Get the derivatives of L
    dLs = jumpoperators_derivatives(Ls_par, theta, dtheta)
   # 2.1 Setup
    njumps = size(traj)[1]
    dpsis = zeros(ComplexF64, sys.NLEVELS, njumps)
    # 2.2 Set up the first jump
    click = traj[1]
    label = click.label
    tau = click.time
    writederivative!(fixlastindex(dpsis, 1),
                     sys.Ls[label], fixlastindex(dLs, label),
                     exp(-1im*tau*sys.Heff),
                     expheff_derivative(Heff_par, tau, theta, dtheta), psi0)
    # In case there are no more jumps, return
    if njumps  == 1
        return dpsis
    end
    # 3. Go over the rest of the jumps
    psitildes = states_atjumps(traj, sys, psi0; normalize=false)
    for k in 2:njumps
      click = traj[k]
      label = click.label
      tau = click.time
      # Calculate the derivative
      writederivative!(fixlastindex(dpsis, k),
                         sys.Ls[label], fixlastindex(dLs, label),
                         exp(-1im*tau*sys.Heff),
                         expheff_derivative(Heff_par, tau, theta, dtheta),
                         fixlastindex(psitildes, k-1), fixlastindex(dpsis, k-1))
   end
   return dpsis

end

function writexi!(xi::SubArray{ComplexF64, 2}, dV::Matrix{ComplexF64},
                  psi::SubArray{ComplexF64, 1}, psi0::Vector{ComplexF64})
    xi .= ((dV*psi0).*adjoint(psi)+psi.*adjoint(dV*psi0))/dot(psi, psi)
end

function writexi!(xi::SubArray{ComplexF64, 2}, V::Matrix{ComplexF64}, dV::Matrix{ComplexF64},
                   psijump::SubArray{ComplexF64, 1}, dpsijump::SubArray{ComplexF64, 1},
                  psi::SubArray{ComplexF64, 1}
                  )
    xi .= ((dV*psijump + V*dpsijump).*adjoint(psi)+psi.*adjoint(dV*psijump + V*dpsijump))/dot(psi, psi)
end

     # tmp = (expheff_derivative(Heff_par, t_given[counter]-t_, theta, dtheta) * psijumps[:, counter_c-1]+
                   # exp(-1im*(t_given[counter]-t_)*sys.Heff)*dpsijumps[:, counter_c-1] )
            # xis[ :, :,counter] = (adjoint(tmp) .* psi[ :, counter] +  adjoint(psi[:, counter]) .* tmp)/dot(psi[:, counter], psi[:, counter])




function monitoringoperator(t_given::Vector{Float64},
    sys::System, Heff_par::Function, Ls_par, traj::Trajectory, psi0::Vector{ComplexF64}, theta::Vector{Float64},
                            dtheta::Vector{Float64})

    # Special case: if the time array is empty, return an empty array
    if isempty(t_given)
        return Array{ComplexF64}(undef, 0, 0, 0)
    end
    psi = states_att(t_given, traj, sys, psi0; normalize=false)
    ntimes = size(t_given)[1]
    njumps = size(traj)[1]
    t_ = 0
    counter = 1
    counter_c = 1
    xis = Array{ComplexF64}(undef, sys.NLEVELS, sys.NLEVELS, ntimes)
    # Edge case
    if isempty(traj)
        while counter <= ntimes
            writexi!(fixlastindex(xis, counter),
                     expheff_derivative(Heff_par, t_given[counter], theta, dtheta),
                     fixlastindex(psi, counter), psi0)
            counter = counter + 1
            if counter > ntimes
                break
            end
        end
        return states
    end
    # Evaluations before first jump
    while (t_given[counter] < traj[counter_c].time) && (counter <= ntimes)
             writexi!(fixlastindex(xis, counter),
                     expheff_derivative(Heff_par, t_given[counter], theta, dtheta),
                     fixlastindex(psi, counter), psi0)
            counter = counter + 1
            if counter > ntimes
                break
            end
    end
    dpsijumps = derivatives_atjumps(sys, Heff_par, Ls_par, traj, psi0, theta, dtheta)
    psijumps = states_atjumps(traj, sys, psi0; normalize=false)
    t_ = t_ + traj[counter_c].time
    counter_c = counter_c + 1
    # Evaluation after first jump
    while (counter_c <= njumps) && (counter <= ntimes)
        timeclick = traj[counter_c].time
        while (t_ < t_given[counter] < t_ + timeclick) && (counter <= ntimes)
            writexi!(fixlastindex(xis, counter), exp(-1im*(t_given[counter]-t_)*sys.Heff),
                     expheff_derivative(Heff_par, t_given[counter]-t_, theta, dtheta),
                     fixlastindex(psijumps, counter_c-1), fixlastindex(dpsijumps, counter_c-1),
                     fixlastindex(psi, counter))
            counter = counter + 1
            if counter > ntimes
                break
            end
         end
       t_ = t_ + timeclick
       counter_c = counter_c + 1
    end

    while counter <= ntimes
        writexi!(fixlastindex(xis, counter), exp(-1im*(t_given[counter]-t_)*sys.Heff),
                     expheff_derivative(Heff_par, t_given[counter]-t_, theta, dtheta),
                     fixlastindex(psijumps, njumps), fixlastindex(dpsijumps, njumps),
                     fixlastindex(psi, counter))

        counter = counter + 1
    end
    return xis
end
