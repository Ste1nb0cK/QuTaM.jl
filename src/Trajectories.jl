module QuTaM

# Dependencies
using LinearAlgebra
import Random
import StatsBase

# Source files
include("structs.jl")
include("functions.jl")

# @reexport using LinearAlgebra  # Re-export all symbols from LinearAlgebra

export run_single_trajectory, run_trajectories, precompute!, SimulParameters, System

end
