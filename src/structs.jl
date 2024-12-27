################# SYSTEM #######################################################
# Note: this struct is make mutable so that it is passed by sharing to the functions
mutable struct System
    NLEVELS::Int64 # Number of levels of the system
    NCHANNELS::Int64 # Number of jump channels
    H::Matrix{ComplexF64} # Hamiltonian
    Ls::Vector{Matrix{ComplexF64}} # List of jump operators
    J::Matrix{ComplexF64} # Sum of Jump operators
    Heff::Matrix{ComplexF64} # Effective Hamiltonian
    # ExpHeff::Function # Function that calculates the exponential of Heff at tau
    # JumpSuperOperators::Vector{Function}  # Jump super operators
    # NoJumpSuperOperator::Function # No-Jump superoperator
    # Inner Constructor
    function System(H::Matrix{ComplexF64}, Ls::Vector{Matrix{ComplexF64}})
        NLEVELS = size(H)[1]
        NCHANNELS = size(Ls)[1] # Number of jump channels
        J = zeros(ComplexF64, NLEVELS, NLEVELS)
        for L in Ls
            J = J + adjoint(L)*L
        end
        CurvyLs = Vector{Function}(undef, NCHANNELS)
       # for k in 1:NCHANNELS
       #     CurvyLs[k] = rho::Matrix{ComplexF64} -> Ls[k]*rho*adjoint(Ls[k])
       # end
        He = H - 0.5im*J
       # expHe(tau::Float64) =  expm(-1im*tau*He)
       # function expcurvyL0(rho::Matrix{ComplexF64}, tau::Float64)
       #     A = expHe(tau)
       #     return A*rho*adjoint(A)
       # end
       # new(NLEVELS, H, Ls, J, He, expHe, CurvyLs, expcurvyL0)
        new(NLEVELS, NCHANNELS, H, Ls, J, He)
    end
end
Base.show(io::IO, s::System) = print(io,
    "System(NLEVELS=$(s.NLEVELS)\nNCHANNELS=$(s.NCHANNELS)\nH=$(s.H)\nLs=$(s.Ls)\nJ=$(s.J))\nHeff=$(s.Heff))")


################ Data Point ################
struct DetectionClick
    time::Float64
    label::Int64
end
const Trajectory = Vector{DetectionClick}
################# SIMULATION PARAMETERS ########################################
# Simulation Struct, it contains the data necessary to run the simulation
# and the precomputed values
mutable struct SimulParameters
    psi0::Vector{ComplexF64}
    nsamples::Int64 # Number of samples in the finegrid
    seed::Int64 # seed
    ntraj::Int64 # Number of trajectories
    multiplier::Float64 # Multiplier to use in the fine grid
    tf::Float64 # Final time
    dt::Float64 # time step for the finegrid
    eps::Float64 # Tolerance for passing WTD normalziation
    function SimulParameters(psi0::Vector{ComplexF64}, tf::Float64,
        s::Int64, ntraj::Int64, nsamples::Int64=10000, m::Float64=10.0,
                             eps::Float64=1e-3)
        deltat = m*tf/nsamples
        new(psi0, nsamples, s, ntraj, m, tf, deltat, eps)
    end
end
Base.show(io::IO, s::SimulParameters) = print(io,
"SimulParameters(psi0=$(s.psi0)\nnsamples=$(s.nsamples)\nseed=$(s.seed)\nntraj=$(s.ntraj))\nmultiplier=$(s.multiplier)\ntf=$(s.tf)\ndt=$(s.dt)\neps=$(s.eps))")
