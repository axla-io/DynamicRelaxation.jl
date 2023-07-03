# A clamped node
function clamped(n_dof)
    @SVector(ones(Bool, n_dof))
end

# node that is completely free
function free(n_dof)
    @SVector(zeros(Bool, n_dof))
end

# a Pinned node that has free rotations but fixed translations
function pinned(n_dof)
    if n_dof == 3
        return clamped(n_dof)
    else
        return SVector{n_dof,Bool}([1, 1, 1, zeros(n_dof - 3)...])
    end
end

const dirmap = Dict(:x => [0, 1, 1], :y => [1, 0, 1], :z => [1, 1, 0])

# a roller node that moves along a cartesian coordinate line and has free rotations
function roller(n_dof, dir)
    SVector{n_dof,Bool}([dirmap[dir]..., zeros(n_dof-3)...])
end

# a roller node that moves on a plane
const planemap = Dict(:xy => [0, 0, 1], :yz => [1, 0, 0], :xz => [0, 1, 0])
function roller_plane(n_dof, dir)
    SVector{n_dof,Bool}([planemap[dir]..., zeros(n_dof-3)...])
end