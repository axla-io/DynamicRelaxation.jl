

function rod_acceleration!(a, x, system::StructuralGraphSystem{Node3DOF}, vertex, s)
    graph = system.graph
    e_map = system.edgemap
    eps = system.elem_props
    x_vert = @view x[:, vertex]
    i_v = UInt8(vertex)
    for neighbor in neighbors(graph, i_v)
        ep = eps[edge_index((i_v, neighbor), e_map)]
        rod_accelerate!(a, x_vert, @view(x[:, neighbor]), ep, s)
    end

    
    return nothing
end


function rod_acceleration!(a, x, system::StructuralGraphSystem{Node6DOF}, vertex, s)
    graph = system.graph
    e_map = system.edgemap
    eps = system.elem_props
    x_vert = @view x[7*(vertex-1)+1:7*vertex]
    i_v = UInt8(vertex)
    for neighbor in neighbors(graph, i_v)
        ep = eps[edge_index((i_v, neighbor), e_map)]
        rod_accelerate!(a, x_vert, @view(x[7*(neighbor-1)+1:7*neighbor]), ep, s, system)
    end

    
    return nothing
end

function f_acceleration!(a, ext_f, i)
    a .+= ext_f[i]
    return nothing
end

function s_min!(s)
    _one = one(eltype(s))
    for i in axes(s, 1)
        s[i] = max(s[i], _one)
    end
    return nothing
end



function rod_accelerate!(a, x0, x1, ep, s, system::StructuralGraphSystem{Node3DOF})
    # Get element length
    element_vec = SVector{3,eltype(x0)}(x1[1] - x0[1], x1[2] - x0[2], x1[3] - x0[3])
    current_length = norm(element_vec)
    rest_length = ep.l_init

    # +++ AXIAL +++
    extension = current_length - rest_length # Unit: [m]

    # +++ FORCES +++
    # Element internal forces
    axial_stiffness = (ep.E * ep.A) / rest_length
    N = axial_stiffness * extension  # Unit: [N]
    a .+= N * element_vec
    s .+= axial_stiffness * abs.(element_vec)

    return nothing

end


function rod_accelerate!(a, x0, x1, ep, s, system::StructuralGraphSystem{Node6DOF})
    # Get element length
    element_vec = SVector{3,eltype(x0)}(x1[1] - x0[1], x1[2] - x0[2], x1[3] - x0[3])
    current_length = norm(element_vec)
    rest_length = ep.l_init

    # +++ AXIAL +++
    extension = current_length - rest_length # Unit: [m]

    # +++ FORCES +++
    # Element internal forces
    axial_stiffness = (ep.E * ep.A) / rest_length
    N = axial_stiffness * extension  # Unit: [N]
    @views a[1:3] .+= N * element_vec
    @views s[1:3] .+= axial_stiffness * abs.(element_vec)

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