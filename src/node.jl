# Struct for storing nodal properties
struct Node3DOF <: NBodySimulator.Body
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
struct CoordinateSystem
    x::SVector{3,Float64}
    y::SVector{3,Float64}
    z::SVector{3,Float64}
end

function CoordinateSystem()
    CoordinateSystem(SVector{3,Float64}(1, 0, 0), SVector{3,Float64}(0, 1, 0), SVector{3,Float64}(0, 0, 1))
end

struct Node6DOF <: NBodySimulator.Body
    # Initial conditions
    r::SVector{3,Float64}
    q::SVector{4,Float64}
    cs::CoordinateSystem
    v::SVector{3,Float64}
    ω::SVector{3,Float64}

    # For constraints
    constrained::Bool
    constraints::SVector{7,Bool}
end

function Node6DOF(pos, constrained, constraints)
    Node6DOF(pos, SVector{4,Float64}(1.0, 0.0, 0.0, 0.0), CoordinateSystem{Float64}(), @SVector(zeros(Float64, 3)), @SVector(zeros(Float64, 3)), constrained, constraints)
end
