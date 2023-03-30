abstract type AbstractGraphSystem <: NBodySimulator.NBodySystem end

struct StructuralGraphSystem{bType<:NBodySimulator.Body} <: NBodySimulator.NBodySystem
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