# Struct for storing nodal properties

abstract type AbstractNode end
struct Node3DOF <: AbstractNode
    # Initial conditions
    r::SVector{3,Float64}
    v::SVector{3,Float64}

    # For constraints
    constrained::Bool
    constraints::SVector{3,Bool}
end

function Node3DOF(pos, constrained, constraints)
    Node3DOF(pos, @SVector(zeros(Float64, 3)), constrained, constraints)
end

function Node3DOF(pos, q0, v0, ω0, constrained, constraints)
    Node3DOF(pos, @SVector(zeros(Float64, 3)), constrained, constraints)
end

num_dof(::Type{Node3DOF}) = 3

struct Node6DOF <: AbstractNode
    # Initial conditions
    r::SVector{3,Float64}
    q::SVector{4,Float64}
    v::SVector{3,Float64}
    ω::SVector{3,Float64}

    # For constraints
    constrained::Bool
    constraints::SVector{7,Bool}
end

function Node6DOF(pos, constrained, constraints)
    Node6DOF(pos, SVector{4,Float64}(1.0, 0.0, 0.0, 0.0), @SVector(zeros(Float64, 3)), @SVector(zeros(Float64, 3)), constrained, constraints)
end

num_dof(::Type{Node6DOF}) = 7
