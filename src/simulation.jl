abstract type StructuralSimulation{T} end

struct RodSimulation{T} <: StructuralSimulation{T}
    system::StructuralGraphSystem{T}
    tspan::Tuple{Float64, Float64}
    dt::Float64
end

function RodSimulation(system::T, tspan, dt) where {T}
    return RodSimulation{eltype(system.bodies)}(system::T, tspan, dt)
end

function get_u0(simulation::T) where T <: StructuralSimulation{Node3DOF}
    system = simulation.system
    bodies = system.bodies
    len = n = length(bodies)
    u_len = 3 * len
    v_len = 3 * len

    u0 = zeros(u_len)
    v0 = zeros(v_len)

    for i in 1:n
        id = 3 * (i - 1) + 1
        u0[id:(id + 2)] = bodies[i].r
        v0[id:(id + 2)] = bodies[i].v
    end

    (u0, v0, n, u_len, v_len)
end

function get_u0(simulation::T) where T <: StructuralSimulation{Node6DOF}
    system = simulation.system
    bodies = system.bodies
    len = n = length(bodies)
    u_len = 7 * len
    v_len = 6 * len

    u0 = zeros(u_len)
    v0 = zeros(v_len)

    for i in 1:n
        bx_id = 7 * (i - 1) + 1
        bv_id = 6 * (i - 1) + 1
        u0[bx_id:(bx_id + 2)] = bodies[i].r
        u0[(bx_id + 3):(bx_id + 6)] = bodies[i].q
        v0[bv_id:(bv_id + 2)] = bodies[i].v
        v0[(bv_id + 3):(bv_id + 5)] = bodies[i].ω
    end

    (u0, v0, n, u_len, v_len)
end

function DiffEqBase.ODEProblem(simulation::S, ext_f) where {T, S<:StructuralSimulation{T}}
    bodies = simulation.system.bodies
    system = simulation.system
    dt = simulation.dt

    (u0, v0, n, u_len, v_len) = get_u0(simulation)
    uv0 = vcat(u0, v0)
    (dx_ids, dr_ids, v_ids, ω_ids) = get_vel_ids(u_len, v_len, system)

    function ode_system!(du, u, p, t)
        # Get positions and element types
        u_v = @view u[1:u_len]
        u_t = eltype(u)
        _z = zero(u_t)

        # Set translation vels
        du[dx_ids] = @view u[v_ids]

        # Set rotation vels
        dr = @view du[dr_ids]
        ω = @view u[ω_ids]

        # Initialize velocity and acceleration
        a = @MVector zeros(u_t, 3)
        s = @MVector zeros(u_t, 3)
        τ = @MVector zeros(u_t, 3)
        j = @MVector zeros(u_t, 3)
        dω = @MVector zeros(u_t, 3)

        @inbounds for i in 1:n
            # Reset accelerations
            reset_accelerations!(a, s, τ, j, dω, _z, simulation)
            body = bodies[i]

            # Accelerate system
            accelerate_system!(a, τ, dω, u_v, system, simulation, body, ext_f, dr, ω, i, s, j, dt, u_t, p)

            # Update accelerations
            update_accelerations!(du, a, dω, u_len, i, simulation)
        end
    end

    return ODEProblem(ode_system!, uv0, simulation.tspan)
end

function reset_accelerations!(a, s, τ, j, dω, _z, simulation::T) where T <: StructuralSimulation{Node6DOF}
    @views a .*= _z
    @views s .*= _z
    @views τ .*= _z
    @views j .*= _z
    @views dω .*= _z
    return nothing
end

function reset_accelerations!(a, s, τ, j, dω, _z, simulation::T) where T <: StructuralSimulation{Node3DOF}
    @views a .*= _z
    @views s .*= _z
    return nothing
end

function update_accelerations!(du, a, dω, u_len, i, simulation::T) where T <: StructuralSimulation{Node6DOF}
    d_id = (u_len) + 6 * (i - 1) + 1
    @views du[d_id:(d_id + 2)] .= a
    @views du[(d_id + 3):(d_id + 5)] .= dω
    return nothing
end

function update_accelerations!(du, a, dω, u_len, i, simulation::T) where T <: StructuralSimulation{Node3DOF}
    d_id = (u_len) + 3 * (i - 1) + 1
    @views du[d_id:(d_id + 2)] .= a
    return nothing
end

function generate_range(n, t1, t2)
    step = (t2 - t1) / (n - 1)
    return [round(Int, t1 + i * step) for i in 0:(n - 1)]
end

function apply_jns!(a, s, dt)
    s .*= dt^2.0 / 2.0
    s_min!(s)
    a .= a ./ s
    return nothing
end

function update_dω!(i, ω, τ, dr, j, dω, u_t, dt)
    dω_id = 3 * (i - 1) + 1
    ω_i = SVector{3, u_t}(ω[dω_id], ω[dω_id] + 1, ω[dω_id] + 2)
    set_rotation_vels!(dr, ω_i, i)

    # Apply moment of inertia
    j .*= dt^2.0 / 2.0
    s_min!(j)
    @views dω .= (τ - scross(ω_i, j .* ω_i)) ./ j
    return nothing
end

function get_vel_ids(u_len, v_len, system::StructuralGraphSystem{Node3DOF})
    dx_ids = get_ids(1, 3, 3, u_len)
    _dr_ids = dx_ids # Dummy variable
    v_ids = get_ids(u_len + 1, 3, 3, u_len + v_len)
    _ω_ids = v_ids # Dummy variable

    return dx_ids, _dr_ids, v_ids, _ω_ids
end

function get_vel_ids(u_len, v_len, system::StructuralGraphSystem{Node6DOF})
    dx_ids = get_ids(1, 3, 7, u_len)
    dr_ids = get_ids(4, 4, 7, u_len)
    v_ids = get_ids(u_len + 1, 3, 6, u_len + v_len)
    ω_ids = get_ids(u_len + 3, 3, 6, u_len + v_len)

    return dx_ids, dr_ids, v_ids, ω_ids
end

function get_state(u, u_len, simulation::T) where T <: StructuralSimulation{Node3DOF}
    x_ids = 1:3:u_len
    y_ids = 2:3:u_len
    z_ids = 3:3:u_len

    state = hcat(u[x_ids], u[y_ids], u[z_ids])'

    return state
end

function get_state(u, u_len, simulation::T) where T <: StructuralSimulation{Node6DOF}
    x_ids = 1:7:u_len
    y_ids = 2:7:u_len
    z_ids = 3:7:u_len

    state = hcat(u[x_ids], u[y_ids], u[z_ids])'

    return state
end

function get_ids(start, step_inc, offset, finish)
    rangelist = Vector{Vector{Int64}}()

    for i in 1:step_inc
        push!(rangelist, collect((start + i - 1):offset:finish))
    end
    return hcat(rangelist...)'[:]
end

function accelerate_system!(a, τ, dω, u_v, system::StructuralGraphSystem{Node6DOF}, simulation::RodSimulation{Node6DOF}, body,
                            ext_f, dr, ω, i, s, j, dt, u_t, p)
    rod_acceleration!(a, τ, u_v, system, body, i, s, j)
    f_acceleration!(a, τ, ext_f, i)
    constrain_acceleration!(a, τ, body)
    apply_jns!(a, s, dt)
    update_dω!(i, ω, τ, dr, j, dω, u_t, dt)
    return nothing
end

function accelerate_system!(a, τ, dω, u_v, system::StructuralGraphSystem{Node3DOF}, simulation::RodSimulation{Node3DOF}, body,
                            ext_f, dr, ω, i, s, j, dt, u_t, p)
    rod_acceleration!(a, u_v, system, i, s)
    f_acceleration!(a, ext_f, i)
    constrain_acceleration!(a, body)
    apply_jns!(a, s, dt)
    return nothing
end
