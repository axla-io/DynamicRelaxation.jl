struct BigonTorqueCondition
    ids::Tuple{Int64, Int64} # Ordered tuple of indices
    cs1::CoordinateSystem # Coordinate system for rod 1
    cs2::CoordinateSystem # Coordinate system for rod 2
    rest_rotation_matrix::SMatrix{3, 3, Float64, 9} # Matrix defining target angle
end

function get_cond_data(cond::BigonTorqueCondition)
    return cond.cs1, cond.cs2, cond.rest_rotation_matrix
end

struct BigonRodSimulation <: StructuralSimulation{Node6DOF}
    system::StructuralGraphSystem{Node6DOF}
    tspan::Tuple{Float64, Float64}
    dt::Float64
    τ_end::Tuple{Int64, Int64}
    τ_start::Tuple{Int64, Int64}
    k::Float64
    nu::Float64
    k_τ::Float64
    nu_τ::Float64
    conds::Vector{BigonTorqueCondition}
end

function get_cs(a, b, e_map, eps)
    ep = eps[edge_index((UInt8(a), UInt8(b)), e_map)]
    return ep.cs
end

function nl_damping_tau(a, b) # A is force and B is damping force
    if sign(a == sign(b))
        return abs(a) > abs(b) ? a - b : a
    else
        return a # We don't want damping to increase force
    end
end

function bigon_forces(a, u, i1, i2, du, simulation)
    x1 = @view u[(7 * (i1 - 1) + 1):(7 * (i1 - 1) + 3)]
    x2 = @view u[(7 * (i2 - 1) + 1):(7 * (i2 - 1) + 3)]

    v1 = @view du[(7 * (i1 - 1) + 1):(7 * (i1 - 1) + 3)]
    v2 = @view du[(7 * (i2 - 1) + 1):(7 * (i2 - 1) + 3)]

    # Calculate contact forces and add to acceleration
    end_distance_vector = x2 - x1
    elastic_force = simulation.k * end_distance_vector
    relative_velocity = v2 - v1
    damping_force = simulation.nu * relative_velocity
    contact_force = elastic_force + damping_force
    a = a + contact_force

    return a
end

function accelerate_system(u_v, system::StructuralGraphSystem{Node6DOF},
                           simulation::BigonRodSimulation, body,
                           ext_f, du, dr_ids, ω, i, dt, u_t, p, t)
    (a, τ, s, j) = rod_acceleration(u_v, system, body, i)

    (a, τ) = f_acceleration(a, τ, ext_f, i, p)
    
    # Apply bigon forces
    if i == simulation.τ_end[1] || i == simulation.τ_end[2]
        i2 = i == simulation.τ_end[1] ? simulation.τ_end[2] : simulation.τ_end[1]
        a = bigon_forces(a, u_v, i, i2, du, simulation)

        cond_id = 1
        τ = τ + bigon_torques(τ, u_v, ω, i, simulation, simulation.conds[cond_id])
        
    elseif i == simulation.τ_start[1] || i == simulation.τ_start[2]
        #i2 = i == simulation.τ_start[1] ? simulation.τ_start[2] : simulation.τ_start[1] # Change  for open behavior
        #a = bigon_forces(a, u_v, i, i2, du, simulation)
        
        cond_id = 2 
        τ = τ + bigon_torques(τ, u_v, ω, i, simulation, simulation.conds[cond_id])
    end

    a = apply_jns!(a, s, dt)
    (a, τ) = constrain_acceleration(a, τ, body)
    dω = update_dω(i, ω, τ, u_v, du, dr_ids, j, u_t, dt)
    return a, dω
end

function bigon_torques(τ, u, ω, i, simulation, cond::BigonTorqueCondition)
    i1, i2 = cond.ids
    # Get coordinate systems and rest rotation matrix
    (cs_1, cs_2, rest_rot_mat) = get_cond_data(cond)

    # Get q1 and q2
    q1 = @view u[(7 * (i1 - 1) + 4):(7 * i1)]
    q2 = @view u[(7 * (i2 - 1) + 4):(7 * i2)]

    # Get ω1 and ω2
    ω1 = @view ω[(3 * (i1 - 1) + 1):(3 * (i1 - 1) + 3)]
    ω2 = @view ω[(3 * (i2 - 1) + 1):(3 * (i2 - 1) + 3)]

    alpha = 1.0
    if i == i1
        alpha = -1.0
    end

    return alpha * apply_bigon_torque(q1, q2, ω1, ω2, rest_rot_mat, cs_1, cs_2, simulation)
end

function apply_bigon_torque(q1, q2, ω1, ω2, rest_rot_mat, cs_1, cs_2, simulation)

    # +++ ROTATIONS +++
    qs_1 = q1[1]
    qs_2 = q2[1]

    qv_1 = SVector{3, eltype(q1)}(q1[2], q1[3], q1[4])
    qv_2 = SVector{3, eltype(q2)}(q2[2], q2[3], q2[4])
    # rel_rot: C_12 = C_1I @ C_I2
    # C_12 is relative rotation matrix from system 1 to system 2
    # C_1I is the rotation from system 1 to the inertial frame (i.e. the world frame)
    # C_I2 is the rotation from the inertial frame to system 2 frame (inverse of system_two_director)
    # rel_rot = system_one_director @ system_two_director.T

    # Get local endplane orientations
    sys1_d = q_vec_rot(qs_1, qv_1, cs_1)
    sys2_d = q_vec_rot(qs_2, qv_2, cs_2)

    rel_rot = sys1_d * sys2_d'

    # error_rot: C_22* = C_21 @ C_12*
    # C_22* is rotation matrix from current orientation of system 2 to desired orientation of system 2
    # C_21 is the inverse of C_12, which describes the relative (current) rotation from system 1 to system 2
    # C_12* is the desired rotation between systems one and two, which is saved in the static_rotation attribute
    dev_rot = rel_rot' * rest_rot_mat

    # compute rotation vectors based on C_22*
    #
    # implementation using custom _inv_rotate
    # rotation vector between identity matrix and C_22*
    rot_vec = inv_rotate(dev_rot')
    #rot_vec[1] *= 1e-3
    #rot_vec[2] = 0.0
    #rot_vec[3] = 0.0

    # rotate rotation vector into inertial frame
    #rot_vec_inertial_frame = sys2_d' * rot_vec
    rot_vec_inertial_frame = sys2_d' * rot_vec

    # deviation in rotation velocity between system 1 and system 2
    # first convert to inertial frame, then take differences
    #dev_omega = sys2_d' * ω2 - sys1_d' * ω1
    dev_omega = ω2 - ω1

    # we compute the constraining torque using a rotational spring - damper system in the inertial frame

    # The opposite torques will be applied to system one and two after rotating the torques into the local frame
    # system_one.external_torques[..., index_one] -= system_one_director @ torque
    # system_two.external_torques[..., index_two] += system_two_director @ torque

    #return simulation.k_τ .* rot_vec_inertial_frame - simulation.nu_τ .* dev_omega
    torque = nl_damping_tau.(simulation.k_τ .* rot_vec_inertial_frame,
                             simulation.nu_τ .* dev_omega) # Use this one for things that work

    return torque
end

function inv_rotate(dev_rotT)
    inv_rot_vec = @MVector zeros(eltype(dev_rotT), 3)

    # Q_{i+i}Q^T_{i} collection
    inv_rot_vec[1] = dev_rotT[3, 2] -
                     dev_rotT[2, 3]

    inv_rot_vec[2] = dev_rotT[1, 3] -
                     dev_rotT[3, 1]

    inv_rot_vec[3] = dev_rotT[2, 1] -
                     dev_rotT[1, 2]

    trace = (dev_rotT[1, 1] +
             dev_rotT[2, 2] +
             dev_rotT[3, 3])

    # TODO HARDCODED bugfix has to be changed. Remove 1e-14 tolerance
    temp = 0.5 * trace - 0.5 - 1e-10
    temp = min(max(temp, -1 + 1e-10), 1 - 1e-10)
    #println(temp)
    theta = acos(temp)

    inv_rot_vec[1] *= -0.5 * theta / sin(theta + 1e-14)
    inv_rot_vec[2] *= -0.5 * theta / sin(theta + 1e-14)
    inv_rot_vec[3] *= -0.5 * theta / sin(theta + 1e-14)

    return inv_rot_vec
end
