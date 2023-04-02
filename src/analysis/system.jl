abstract type AbstractGraphSystem <: NBodySimulator.NBodySystem end

struct StructuralGraphSystem{bType<:NBodySimulator.Body} <: AbstractGraphSystem
    # nbody
    bodies::Vector{bType}

    # graph with same order as nbody
    graph::StaticGraph{UInt8,UInt8} #TODO: best way to define graph types?

    # custom
    elem_props::Vector{ElementProperties{Float64}} # Also check how to properly set element properties
    edgemap::Dict{Tuple{UInt8,UInt8},Int}
end

function edge_index(e::Tuple{UInt8,UInt8}, e_map::Dict{Tuple{UInt8,UInt8},Int})
    res = 0
    rev_e = reverse(e)
    if haskey(e_map, e)
        res = e_map[e]
    elseif haskey(e_map, rev_e)
        res = e_map[rev_e]
    end
    return res
end



function default_system(graph, node_type, sys_type)

    if node_type == Node3DOF
        n_dof = 3 

    elseif node_type == Node6DOF
        n_dof = 7
    else
        error("invalid node type")
    end

    E = 210 * 1e9               # [Pa]
    Iy = Iz = 4.7619 * 1e-7     # [m^4]
    A = 4.7619 * 1e-4           # [m^2]
    G = 78 * 1e9                # [Pa]
    It = 2 * Iy                 # [m^4]
    l_init = 1.0

    ep = ElementProperties{Float64}(E, A, Iy, Iz, G, It, l_init)

    if sys_type == :catenary # Is equivalent to simply supported beam for beam structure
        #@infiltrate
        np_fix = node_type{Float64}(@SVector(zeros(3)), true, pinned(n_dof))
        np_free = [node_type{Float64}(SVector{3,Float64}([i - 1, 0.0, 0.0]), false, free(n_dof)) for i in 2:n_pt-1]
        nodes = vcat(np_fix, np_free..., node_type{Float64}(SVector{3,Float64}([n_elem, 0.0, 0.0]), true, pinned(n_dof))) # Assuming same order as in graph

    elseif sys_type == :cantilever
        np_fix = node_type{Float64}(@SVector(zeros(3)), true, clamped(n_dof))
        np_free = [node_type{Float64}(SVector{3,Float64}([i - 1, 0.0, 0.0]), false, free(n_dof)) for i in 2:n_pt]
        nodes = vcat(np_fix, np_free...) # Assuming same order as in graph

    elseif sys_type == :elastica
        np_fix = node_type{Float64}(@SVector(zeros(3)), true, pinned(n_dof))
        np_free = [node_type{Float64}(SVector{3,Float64}([i - 1, 0.0, 0.0]), false, free(n_dof)) for i in 2:n_pt-1]
        np_roller = node_type{Float64}(SVector{3,Float64}([n_elem, 0.0, 0.0]), true, roller(n_dof, :x))
        nodes = vcat(np_fix, np_free..., np_roller) # Assuming same order as in graph

    else
        error("invalid system type")
    end
    eps = [ep for _e in edges(graph)] # Assuming same order as in graph

    edgelist = collect(edges(graph))

    edgemap = Dict{Tuple{UInt8,UInt8},Int}((src(e), dst(e)) => i for (i, e) in enumerate(edgelist))

    StructuralGraphSystem{node_type}(nodes, graph, eps, edgemap)

end