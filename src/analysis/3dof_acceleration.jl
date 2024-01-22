function rod_acceleration(x, system::StructuralGraphSystem{Node3DOF}, vertex,)
    graph = system.graph
    e_map = system.edgemap
    eps = system.elem_props
    v_i = 3 * (vertex - 1) + 1
    x_vert = @view x[v_i:(v_i + 2)]
    #i_v = UInt8(vertex)
    i_v = vertex
    a = @SVector zeros(eltype(x), 3)
    s = @SVector zeros(eltype(x), 3)

    for neighbor in neighbors(graph, i_v)
        n_i = 3 * (neighbor - 1) + 1
        ep = eps[edge_index((i_v, neighbor), e_map)]
        (a, s) = rod_accelerate(a, x_vert, @view(x[n_i:(n_i + 2)]), ep, s)
    end
    return a, s
end

function f_acceleration(a, ext_f, i::Int)
    a = SA[a[1] + ext_f[i][1], a[2] + ext_f[i][2], a[3] + ext_f[i][3]]
    return a
end

function s_min!(s)
    _one = one(eltype(s))
    return SA[max(s[1], _one), max(s[2], _one), max(s[3], _one)]
end

function rod_accelerate(a, x0, x1, ep, s)
    # Get element length
    element_vec = SVector{3, eltype(x0)}(x1[1] - x0[1], x1[2] - x0[2], x1[3] - x0[3])
    current_length = norm(element_vec)
    rest_length = ep.l_init

    axial_stiffness = (ep.E * ep.A) / rest_length
    # Deformation
    extension = current_length - rest_length # Unit: [m]
    if !(ep.cable) || ep.cable && extension >= -1e-7
        # Element internal forces
        N = axial_stiffness * extension  # Unit: [N]
        a = a + N * element_vec
    end
    s = s + axial_stiffness * abs.(element_vec)
    return a, s
end

function get_constrained(constraints::SVector{3, Bool}, x1::S) where {T, S <: SVector{3, T}}
    x1_c = SVector{3, T}(constraints[i] ? (i < 4 ? zero(T) : x1[i]) : x1[i]
                         for i in axes(x1, 1))
    return x1_c
end

function constrain_acceleration(a, body)
    if body.constrained == true
        constraints = body.constraints

        a = get_constrained(constraints, a)
    end
    return a
end
