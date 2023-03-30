

function rod_acceleration!(a, x, system, vertex)
    graph = system.graph
    e_map = system.edgemap
    eps = system.elem_props
    x_vert = @view x[:, vertex]

    for neighbor in neighbors(graph, vertex)
        ep = @view eps[edge_index((vertex, neighbor), e_map)]
        rod_accelerate!(a, x_vert, @view x[:, neighbor], ep)
    end

    return nothing
end







function f_acceleration!(a, ext_f, i)
    a += ext_f[i]
    return nothing
end

function rod_accelerate!(a, x0, x1, ep)
    # Get element length
    element_vec = SVector{3,eltype(x0)}(x0 .- x1)
    current_length = norm(element_vec)
    rest_length = ep.l_init

    # +++ AXIAL +++
    extension = current_length - rest_length # Unit: [m]

    # +++ FORCES +++
    # Element internal forces
    axial_stiffness = (ep.E * ep.A) / rest_length
    N = axial_stiffness * extension  # Unit: [N]

    a += N * element_vec

    return nothing

end