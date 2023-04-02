using DiffEqCallbacks
using Graphs
using StaticArrays
using LinearAlgebra
using NBodySimulator
using StaticGraphs
using ForwardDiff

using Infiltrator
using JET

include("elem.jl")
include("node.jl")
include("simulation.jl")
include("analysis/loads.jl")
include("analysis/rotations.jl")
include("analysis/constraints.jl")
include("analysis/3dof_acceleration.jl")
include("analysis/6dof_acceleration.jl")
include("analysis/system.jl")
include("optimization/load_finding.jl")