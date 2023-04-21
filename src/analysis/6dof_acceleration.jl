get_properties(ep) = (ep.E, ep.A, ep.Iy, ep.Iz, ep.G, ep.It)

function rod_accelerate!(a1, a2, τ1, τ2, u0, u1, body_i, body_j, ep, s1, s2, j1, j2)
    # Get element properties
    (E, A, Iy, Iz, G, It) = get_properties(ep)

    # Get element length
    element_vec = SVector{3, eltype(u0)}(u1[1] - u0[1], u1[2] - u0[2], u1[3] - u0[3])
    current_length = norm(element_vec)
    rest_length = ep.l_init

    # +++ ROTATIONS +++
    qs_i = u0[4]
    qs_j = u1[4]

    qv_i = SVector{3, eltype(u0)}(u0[5], u0[6], u0[7])
    qv_j = SVector{3, eltype(u0)}(u1[5], u1[6], u1[7])

    cs_i = body_i.cs
    cs_j = body_j.cs

    # Get local endplane orientations
    y0 = q_vec_rot(qs_i, qv_i, cs_i.y)
    z0 = q_vec_rot(qs_i, qv_i, cs_i.z)
    x0 = q_vec_rot(qs_i, qv_i, cs_i.x)

    y1 = q_vec_rot(qs_j, qv_j, cs_j.y)
    z1 = q_vec_rot(qs_j, qv_j, cs_j.z)
    x1 = q_vec_rot(qs_j, qv_j, cs_j.x)

    #Bending angle changes around local axes
    inv_current_length = 1.0 / current_length
    theta_y0 = (z0 ⋅ element_vec) * inv_current_length
    theta_z0 = -(y0 ⋅ element_vec) * inv_current_length #NB! Negative sign

    theta_y1 = (z1 ⋅ element_vec) * inv_current_length
    theta_z1 = -(y1 ⋅ element_vec) * inv_current_length #NB! Negative sign

    #Twist angle change around element axis
    theta_x = ((y0 ⋅ z1) - (y1 ⋅ z0)) * 0.5

    # +++ AXIAL +++
    r_30 = rest_length / 30.0
    ext_a = (current_length^2 - rest_length^2) / (2.0 * rest_length)

    ext_b = (r_30 * 0.5) * (4.0 * (theta_y0^2 + theta_z0^2) -
             2.0 * ((theta_y0 * theta_y1)
                    +
                    (theta_z0 * theta_z1)) + 4.0 * (theta_y1^2 + theta_z1^2))
    extension = ext_a + ext_b # Unit: [m]

    # +++ FORCES +++
    # Element internal forces
    inv_rest_length = 1.0 / rest_length
    axial_stiffness = (E * A) * inv_rest_length
    N = axial_stiffness * extension  # Unit: [N]

    # +++ MOMENTS +++
    M_y0 = ((N * r_30) * ((4.0 * theta_y0) - theta_y1)) +
           (((E * Iy) * inv_rest_length) * ((4.0 * theta_y0) + (2.0 * theta_y1)))           #Unit: [Nm]
    M_z0 = ((N * r_30) * ((4.0 * theta_z0) - theta_z1)) +
           (((E * Iz) * inv_rest_length) * ((4.0 * theta_z0) + (2.0 * theta_z1)))           #Unit: [Nm]

    M_y1 = ((N * r_30) * ((4.0 * theta_y1) - theta_y0)) +
           (((E * Iy) * inv_rest_length) * ((4.0 * theta_y1) + (2.0 * theta_y0)))           #Unit: [Nm]
    M_z1 = ((N * r_30) * ((4.0 * theta_z1) - theta_z0)) +
           (((E * Iz) * inv_rest_length) * ((4.0 * theta_z1) + (2.0 * theta_z0)))           #Unit: [Nm]

    M_x = ((G * It) * inv_rest_length) * theta_x            #Unit: [Nm]

    #Force start
    F0_x = inv_rest_length *
           ((N * element_vec[1]) + (M_y0 * z0[1]) - (M_z0 * y0[1]) + (M_y1 * z1[1]) -
            (M_z1 * y1[1]))
    F0_y = inv_rest_length *
           ((N * element_vec[2]) + (M_y0 * z0[2]) - (M_z0 * y0[2]) + (M_y1 * z1[2]) -
            (M_z1 * y1[2]))
    F0_z = inv_rest_length *
           ((N * element_vec[3]) + (M_y0 * z0[3]) - (M_z0 * y0[3]) + (M_y1 * z1[3]) -
            (M_z1 * y1[3]))

    #Moment start
    #i=1, j=2, k=3
    M0x_pos = -(((M_y0 * element_vec[3] * z0[2]) * inv_rest_length) -
                ((M_z0 * element_vec[3] * y0[2]) * inv_rest_length) +
                ((M_x * ((y0[2] * z1[3]) - (z0[2] * y1[3]))) * 0.5))

    #i=1, j=3, k=2
    M0x_neg = (((M_y0 * element_vec[2] * z0[3]) * inv_rest_length) -
               ((M_z0 * element_vec[2] * y0[3]) * inv_rest_length) +
               ((M_x * ((y0[3] * z1[2]) - (z0[3] * y1[2]))) * 0.5))

    #i=2, j=3, k=1
    M0y_pos = -(((M_y0 * element_vec[1] * z0[3]) * inv_rest_length) -
                ((M_z0 * element_vec[1] * y0[3]) * inv_rest_length) +
                ((M_x * ((y0[3] * z1[1]) - (z0[3] * y1[1]))) * 0.5))

    #i=2, j=1, k=3
    M0y_neg = (((M_y0 * element_vec[3] * z0[1]) * inv_rest_length) -
               ((M_z0 * element_vec[3] * y0[1]) * inv_rest_length) +
               ((M_x * ((y0[1] * z1[3]) - (z0[1] * y1[3]))) * 0.5))

    #i=3, j=1, k=2
    M0z_pos = -(((M_y0 * element_vec[2] * z0[1]) * inv_rest_length) -
                ((M_z0 * element_vec[2] * y0[1]) * inv_rest_length) +
                ((M_x * ((y0[1] * z1[2]) - (z0[1] * y1[2]))) * 0.5))

    #i=3, j=2, k=1
    M0z_neg = (((M_y0 * element_vec[1] * z0[2]) * inv_rest_length) -
               ((M_z0 * element_vec[1] * y0[2]) * inv_rest_length) +
               ((M_x * ((y0[2] * z1[1]) - (z0[2] * y1[1]))) * 0.5))

    #Moment end
    #i=1, j=2, k=3
    M1x_pos = -(((M_y1 * element_vec[3] * z1[2]) * inv_rest_length) -
                ((M_z1 * element_vec[3] * y1[2]) * inv_rest_length) -
                ((M_x * ((y0[2] * z1[3]) - (z0[2] * y1[3]))) * 0.5))

    #i=1, j=3, k=2
    M1x_neg = (((M_y1 * element_vec[2] * z1[3]) * inv_rest_length) -
               ((M_z1 * element_vec[2] * y1[3]) * inv_rest_length) -
               ((M_x * ((y0[3] * z1[2]) - (z0[3] * y1[2]))) * 0.5))

    #i=2, j=3, k=1
    M1y_pos = -(((M_y1 * element_vec[1] * z1[3]) * inv_rest_length) -
                ((M_z1 * element_vec[1] * y1[3]) * inv_rest_length) -
                ((M_x * ((y0[3] * z1[1]) - (z0[3] * y1[1]))) * 0.5))

    #i=2, j=1, k=3
    M1y_neg = (((M_y1 * element_vec[3] * z1[1]) * inv_rest_length) -
               ((M_z1 * element_vec[3] * y1[1]) * inv_rest_length) -
               ((M_x * ((y0[1] * z1[3]) - (z0[1] * y1[3]))) * 0.5))

    #i=3, j=1, k=2
    M1z_pos = -(((M_y1 * element_vec[2] * z1[1]) * inv_rest_length) -
                ((M_z1 * element_vec[2] * y1[1]) * inv_rest_length) -
                ((M_x * ((y0[1] * z1[2]) - (z0[1] * y1[2]))) * 0.5))

    #i=3, j=2, k=1
    M1z_neg = (((M_y1 * element_vec[1] * z1[2]) * inv_rest_length) -
               ((M_z1 * element_vec[1] * y1[2]) * inv_rest_length) -
               ((M_x * ((y0[2] * z1[1]) - (z0[2] * y1[1]))) * 0.5))

    # Update forces
    a1[1] += F0_x
    a1[2] += F0_y
    a1[3] += F0_z

    a2[1] -= F0_x
    a2[2] -= F0_y
    a2[3] -= F0_z

    τ1[1] += M0x_pos + M0x_neg
    τ1[2] += M0y_pos + M0y_neg
    τ1[3] += M0z_pos + M0z_neg

    τ2[1] += M1x_pos + M1x_neg
    τ2[2] += M1y_pos + M1y_neg
    τ2[3] += M1z_pos + M1z_neg

    # Update stiffnesses
    s1 .+= axial_stiffness * abs.(element_vec)
    s2 .+= axial_stiffness * abs.(element_vec)
    update_j!(y0 + y1, z0 + z1, x0 + x1, j1, E, Iy, Iz, G, It, inv_rest_length)
    update_j!(y0 + y1, z0 + z1, x0 + x1, j2, E, Iy, Iz, G, It, inv_rest_length)

    return nothing
end

function update_j!(y_m, z_m, x_m, j, E, Iy, Iz, G, It, inv_rest_length)
    # Update stiffness
    j .+= (E .* (Iy .* abs.(y_m) + Iz .* abs.(z_m)) + G .* It .* abs.(x_m)) *
          inv_rest_length

    return nothing
end

function constrain_acceleration!(a, τ, body)
    if body.constrained == true
        constraints = body.constraints
        _zero = zero(eltype(a))
        for i in 1:length(body.constraints)
            if constraints[i] == true
                if i < 4
                    a[i] = _zero
                else
                    τ[i] = _zero
                end
            end
        end
    end
    return nothing
end

function f_acceleration!(a, τ, ext_f, i)
    for j in 1:3
        a[j] += ext_f[i][j]
        τ[j] += ext_f[i][j + 3]
    end
    return nothing
end

function rod_acceleration!(a1, a2, τ1, τ2, x, system::StructuralGraphSystem{Node6DOF}, ep,
                           id1, id2, body1, body2, s1, s2, j1, j2)
    v_1 = 3 * (id1 - 1) + 1
    x_1 = @view x[v_1:(v_1 + 2)]

    v_2 = 3 * (id2 - 1) + 1
    x_2 = @view x[v_2:(v_2 + 2)]

    rod_accelerate!(a1, a2, τ1, τ2, x_1, x_2, body1, body2, ep, s1, s2, j1, j2)

    return nothing
end
