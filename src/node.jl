# Struct for storing nodal properties

struct NodeProperties{T}
    # For time integration
    M::SMatrix{3, 3, T}
    M_inv::SMatrix{3, 3, T}
    J::SMatrix{3, 3, T}
    J_inv::SMatrix{3, 3, T} 

    # For constraints
    constrained::Bool
    constraints::SVector{7, Bool}
end