function set_rotation_vels!(u::AbstractVector{T1}, du::AbstractVector{T2}, dr_ids, ω_i, i) where {T1 <: Real, T2 <: Real}
    dr_id = 4 * (i - 1) + 1

    #r_i = @view u[@view(dr_ids[dr_id:(dr_id + 3)])]
    # r_i_s = SVector{4, T1}(r_i[1], r_i[2], r_i[3], r_i[4])
    #r_i = @view u[@view(dr_ids[dr_id:(dr_id + 3)])]
    r_i_s = SVector{4, T1}(u[dr_ids[dr_id]], u[dr_ids[dr_id+1]], u[dr_ids[dr_id+2]], u[dr_ids[dr_id+3]])
    
    # Do rotation according to Rucker, “Integrating Rotations Using Nonunit Quaternions.”
    dr_i = @view du[@view(dr_ids[dr_id:(dr_id + 3)])]

    dr_i .= f_q_dot(r_i_s, ω_i)
    return nothing
end

function scross(a::AbstractVector{T1}, b::AbstractVector{T2}) where {T1 <: Real, T2 <: Real}
    a1, a2, a3 = a
    b1, b2, b3 = b

    SVector{3, T1}(a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1)
end

# Ported from https://github.com/eayvali/Integrating-Rigid-Body-Rotations/blob/master/code/f_q_dot.m
function f_q_dot(q::SVector{4, T}, w::SVector{3, T}; k = 0.1) where {T <: Real}
    W = get_W(w) # Skew-symmetric matrix of angular velocities
    c = k * (1.0 .- q' * q)
    return 0.5 * W * q + c .* q # time derivative of quaternion
end

# Ported from https://github.com/eayvali/Integrating-Rigid-Body-Rotations/blob/master/code/f_q_dot.m
function get_W(w::SVector{3, T}) where {T <: Real}
    _z = zero(T)
    W = SA[_z -w[1] -w[2] -w[3]
           w[1] _z w[3] -w[2]
           w[2] -w[3] _z w[1]
           w[3] w[2] -w[1] _z]
    return W
end

# Quaternion rotation of vector, from: https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion

function q_vec_rot(qs, qv, v)
    return 2.0 * dot(qv, v) .* qv +
           (qs^2 - dot(qv, qv)) .* v +
           2.0 * qs .* scross(qv, v)
end
