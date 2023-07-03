module DynamicRelaxation

# Arrays
using StaticArrays
using LinearAlgebra

# Graph deps
using Graphs
using StaticGraphs

# Differential equation solving
using DiffEqBase
using DiffEqCallbacks
using NBodySimulator
using SteadyStateDiffEq
using ForwardDiff

include("elem.jl")
include("node.jl")
include("analysis/system.jl")
include("simulation.jl")
include("analysis/loads.jl")
include("analysis/rotations.jl")
include("analysis/constraints.jl")
include("analysis/callbacks.jl")
include("analysis/3dof_acceleration.jl")
include("analysis/6dof_acceleration.jl")
include("optimization/load_finding.jl")

# Elements and nodes
export ElementProperties, CoordinateSystem, Node3DOF, Node6DOF

# Callbacks
export velocitydecay!, velocityreset!, ke_condition, ke_termination_cond

# Constraints
export clamped, free, pinned, roller

# Loads
export Px, Py, Pz, Mx, My, Mz, uniform_load, point_loads

# System 
export StructuralGraphSystem, default_system

# Simulation
export LoadScaleRodSimulation, RodSimulation, get_u0, get_vel_ids, get_state

# Plotting
export generate_range

end  # module
