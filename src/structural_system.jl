abstract type AbstractGraphSystem <: NBodySimulator.NBodySystem end

struct StructuralGraphSystem{bType<:NBodySimulator.MassBody,fType<:Real} <: NBodySimulator.BasicPotentialSystem
    # nbody
    bodies::Vector{bType}
    ext_f::Vector{SVector{3, fType}}
    
    # graph with same order as nbody
    graph::StaticGraph{UInt8, UInt8} #TODO: best way to define graph types?

    # custom
    node_props::Vector{NodeProperties{Float64}} # Double check if this is actually the same as bodies
    elem_props::Vector{ElementProperties{Float64}} # Also check how to properly set element properties
end
