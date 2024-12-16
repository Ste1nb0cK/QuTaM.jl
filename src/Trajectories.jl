module Trajectories

# Dependencies
using LinearAlgebra
import Random
import StatsBase

# Source files
include("functions.jl")
include("structs.jl")

# @reexport using LinearAlgebra  # Re-export all symbols from LinearAlgebra

export run_single_trajectory, run_trajectories, SimulParameters

end
