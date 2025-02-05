module BackAction

# Dependencies
using LinearAlgebra
using Statistics
using ProgressMeter
using Base.Threads
import Random
import StatsBase

# Source files
include("structs.jl")
include("precompute.jl")
include("functions_jump.jl")
include("run_trajectories.jl")
include("monitoring.jl")
# Utilities
include("../util/pauli_m.jl")
include("../util/rd_ex.jl")
include("../util/rd_temperature_ex.jl")
include("../util/rf_ex.jl")
include("../util/rk4.jl")
export
    # Structs
    System, SimulParameters, DetectionClick, Trajectory,
    # Functions
    run_trajectories,
    sample_single_trajectory,
    states_att,
    states_atjumps,
    monitoringoperator
end
