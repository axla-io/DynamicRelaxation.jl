# Struct for storing nodal properties

struct Node3DOF{T}<:NBodySimulator.Body
    # Initial conditions
    r::SVector{3, T}
    v::SVector{3, T}

    # For time integration
    M::SMatrix{3, 3, T}
    M_inv::SMatrix{3, 3, T}

    # For constraints
    constrained::Bool
    constraints::SVector{3, Bool}
end

struct Node6DOF{T}<:NBodySimulator.Body
    # Initial conditions
    r::SVector{3, T}
    q::SVector{3, T}
    v::SVector{3, T}
    ω::SVector{3, T}

    # For time integration
    M::SMatrix{3, 3, T}
    M_inv::SMatrix{3, 3, T}
    J::SMatrix{3, 3, T}
    J_inv::SMatrix{3, 3, T} 

    # For constraints
    constrained::Bool
    constraints::SVector{7, Bool}
end