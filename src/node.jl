# Struct for storing nodal properties

struct Node3DOF{T}<:NBodySimulator.Body
    # Initial conditions
    r::SVector{3, T}
    v::SVector{3, T}

    # For constraints
    constrained::Bool
    constraints::SVector{3, Bool}
end

function Node3DOF{T}(pos, constrained, constraints) where T <:Real
    Node3DOF{T}(pos, @SVector(zeros(T, 3)), constrained, constraints)
end

struct Node6DOF{T}<:NBodySimulator.Body
    # Initial conditions
    r::SVector{3, T}
    q::SVector{4, T}
    v::SVector{3, T}
    ω::SVector{3, T}

    # For constraints
    constrained::Bool
    constraints::SVector{7, Bool}
end


function Node6DOF{T}(pos, constrained, constraints) where T <:Real
    Node6DOF{T}(pos, SVector{4,T}(1.0, 0.0, 0.0,0.0), @SVector(zeros(T, 3)), @SVector(zeros(T, 3)), constrained, constraints)
end