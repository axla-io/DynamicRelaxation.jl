
Px(c, system::StructuralGraphSystem{Node6DOF}) = c * [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Py(c, system::StructuralGraphSystem{Node6DOF}) = c * [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
Pz(c, system::StructuralGraphSystem{Node6DOF}) = c * [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]

Px(c, system::StructuralGraphSystem{Node3DOF}) = c * [1.0, 0.0, 0.0]
Py(c, system::StructuralGraphSystem{Node3DOF}) = c * [0.0, 1.0, 0.0]
Pz(c, system::StructuralGraphSystem{Node3DOF}) = c * [0.0, 0.0, 1.0]

Mx(c, system::StructuralGraphSystem{Node6DOF}) = c * [ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
My(c, system::StructuralGraphSystem{Node6DOF}) = c * [ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
Mz(c, system::StructuralGraphSystem{Node6DOF}) = c * [ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]


function uniform_load(v::Vector{T}, system) where {T<:Real}
    N = length(v)
    return [SVector{N,T}(v) for i in 1:Int(nv(system.graph))]
end


function point_loads(vs, is, system)
    if length(vs) > 0
        v_dim = length(vs[1])
    else
        error("no loads")
    end

    loads = Vector{SVector{v_dim,Float64}}()
    for i in 1:Int(nv(system.graph))
        if i in is
            #@infiltrate
            load = SVector{v_dim,Float64}(vs[findfirst(isequal(i), is)])
        else
            load = @SVector zeros(v_dim)
        end
        push!(loads, load)
    end
    return loads
end