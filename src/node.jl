# Struct for storing nodal properties
struct Node3DOF{T} <: NBodySimulator.Body
    # Initial conditions
    r::SVector{3,T}
    v::SVector{3,T}

    # For constraints
    constrained::Bool
    constraints::SVector{3,Bool}
end

function Node3DOF{T}(pos, constrained, constraints) where {T<:Real}
    Node3DOF{T}(pos, @SVector(zeros(T, 3)), constrained, constraints)
end
struct CoordinateSystem{T}
    x::SVector{3,T}
    y::SVector{3,T}
    z::SVector{3,T}
end

function CoordinateSystem{T}()where T <:Real
    CoordinateSystem{T}(SVector{3,T}(1, 0, 0), SVector{3,T}(0, 1, 0), SVector{3,T}(0, 0, 1))
end

struct Node6DOF{T} <: NBodySimulator.Body
    # Initial conditions
    r::SVector{3,T}
    q::SVector{4,T}
    cs::CoordinateSystem{T}
    v::SVector{3,T}
    ω::SVector{3,T}

    # For constraints
    constrained::Bool
    constraints::SVector{7,Bool}
end

function Node6DOF{T}(pos, constrained, constraints) where {T<:Real}
    Node6DOF{T}(pos, SVector{4,T}(1.0, 0.0, 0.0, 0.0), CoordinateSystem{T}(), @SVector(zeros(T, 3)), @SVector(zeros(T, 3)), constrained, constraints)
end
