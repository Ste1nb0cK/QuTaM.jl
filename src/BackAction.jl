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
export
    # Structs
    System, SimulParameters, DetectionClick, Trajectory,
    # Functions
    run_trajectories,
    sample_single_trajectory,
    evaluate_at_t,
    states_at_jumps,
    MonitoringOperator
end
