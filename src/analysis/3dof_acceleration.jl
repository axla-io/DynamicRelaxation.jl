function rod_acceleration!(a1, a2, x, system::StructuralGraphSystem{Node3DOF}, ep, id1, id2, s1, s2)
    
    v_1 = 3*(id1 - 1) + 1
    x_1 = @view x[v_1:v_1+2]

    v_2 = 3*(id2 - 1) + 1
    x_2 = @view x[v_2:v_2+2]

    rod_accelerate!(a1, a2, x_1, x_2, ep, s1, s2)

    return nothing
end

function f_acceleration!(a, ext_f, i)
    for j = 1:3
        a[j] += ext_f[i][j]
    end
    return nothing
end

function s_min!(s)
    _one = one(eltype(s))
    for i in axes(s, 1)
        s[i] = max(s[i], 0.5_one)
    end
    return nothing
end

function rod_accelerate!(a1, a2, x0, x1, ep, s1, s2)
    # Get element length
    element_vec = SVector{3,eltype(x0)}(x1[1] - x0[1], x1[2] - x0[2], x1[3] - x0[3])
    current_length = norm(element_vec)
    rest_length = ep.l_init

    # Deformation
    extension = current_length - rest_length # Unit: [m]

    # Element internal forces
    axial_stiffness = (ep.E * ep.A) / rest_length
    N = axial_stiffness * extension  # Unit: [N]
    F = N * element_vec
    a1 .+= F
    a2 .-= F

    s1 .+= axial_stiffness * abs.(element_vec)
    s2 .+= axial_stiffness * abs.(element_vec)
    return nothing

end

function constrain_acceleration!(a, body)
    if body.constrained == true
        constraints = body.constraints
        _zero = zero(eltype(a))
        for i = 1:length(body.constraints)
            if constraints[i] == true
                a[i] = _zero
            end
        end
    end
    return nothing
end
