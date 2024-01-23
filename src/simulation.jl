abstract type StructuralSimulation{T} end

struct RodSimulation{T} <: StructuralSimulation{T}
    system::StructuralGraphSystem{T}
    tspan::Tuple{Float64, Float64}
    dt::Float64
end

function RodSimulation(system::T, tspan, dt) where {T}
    return RodSimulation{eltype(system.bodies)}(system::T, tspan, dt)
end

function get_u0(simulation::T) where {T <: StructuralSimulation{Node3DOF}}
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

function get_u0(simulation::T) where {T <: StructuralSimulation{Node6DOF}}
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

interpol(x) = x < 0.1 ? 10x : one(x)

function DiffEqBase.ODEProblem(simulation::S, ext_f; constrain_jac = true) where {T, S <: StructuralSimulation{T}}
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

        # Set translation vels
        @views du[dx_ids] .= u[v_ids]

        # Get angular velocities
        ω = @view u[ω_ids]

        t_mult = 1.0
        
        if typeof(p) == typeof(t_mult) # This is a hack and should be changed later
            t_mult *= p
        end

        @inbounds for i in 1:n

            # Get current body
            body = bodies[i]
            
            # Accelerate system
            (a, dω) = accelerate_system(u_v, system, simulation, body, ext_f, du,
                                        dr_ids, ω, i, dt, u_t, t_mult, t)
            # Update accelerations
            update_accelerations!(du, a, dω, u_len, i, simulation)
        end
    end

    ode_jac = constrain_jac ? get_ode_jac(ode_system!, u_len, simulation)  : nothing
    #ode_jac_prototype = get_jac_prototype(system, u_len, v_len)
    ode_jac_prototype = nothing
    ode_f = ODEFunction(ode_system!, jac = ode_jac, jac_prototype = ode_jac_prototype)

    return ODEProblem(ode_f, uv0, simulation.tspan)
end

function update_accelerations!(du, a, dω, u_len, i,
                               simulation::T) where {T <: StructuralSimulation{Node6DOF}}
    d_id = (u_len) + 6 * (i - 1) + 1
    @views du[d_id:(d_id + 2)] .= a
    @views du[(d_id + 3):(d_id + 5)] .= dω
    return nothing
end

function update_accelerations!(du, a, dω, u_len, i,
                               simulation::T) where {T <: StructuralSimulation{Node3DOF}}
    d_id = (u_len) + 3 * (i - 1) + 1
    @views du[d_id:(d_id + 2)] .= a
    return nothing
end

function generate_range(n, t1, t2)
    step = (t2 - t1) / (n - 1)
    return [round(Int, t1 + i * step) for i in 0:(n - 1)]
end

function apply_jns!(a, s, dt)
    s = s * dt^2.0 / 2.0
    s = s_min!(s)
    a = a ./ s
    return a
end

function update_dω(i, ω, τ, u_v, du, dr_ids, j, u_t, dt)
    dω_id = 3 * (i - 1) + 1
    ω_i = SA[ω[dω_id], ω[dω_id + 1], ω[dω_id + 2]]
    set_rotation_vels!(u_v, du, dr_ids, ω_i, i)

    # Apply moment of inertia
    j = j * dt^2.0 / 2.0
    j = s_min!(j)
    dω = (τ - scross(ω_i, j .* ω_i)) ./ j
    return dω
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
    ω_ids = get_ids(u_len + 4, 3, 6, u_len + v_len)

    return dx_ids, dr_ids, v_ids, ω_ids
end

function get_state(u, u_len, simulation::T) where {T <: StructuralSimulation{Node3DOF}}
    x_ids = 1:3:u_len
    y_ids = 2:3:u_len
    z_ids = 3:3:u_len

    state = hcat(u[x_ids], u[y_ids], u[z_ids])'

    return state
end

function get_state(u, u_len, simulation::T) where {T <: StructuralSimulation{Node6DOF}}
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

function accelerate_system(u_v, system::StructuralGraphSystem{Node6DOF},
                           simulation::RodSimulation{Node6DOF}, body,
                           ext_f, du, dr_ids, ω, i, dt, u_t, p, t)
    (a, τ, s, j) = rod_acceleration(u_v, system, body, i)

    (a, τ) = f_acceleration(a, τ, ext_f, i, p)
    (a, τ) = constrain_acceleration(a, τ, body)
    a = apply_jns!(a, s, dt)
    dω = update_dω(i, ω, τ, u_v, du, dr_ids, j, u_t, dt)
    return a, dω
end

function accelerate_system(u_v, system::StructuralGraphSystem{Node3DOF},
                           simulation::RodSimulation{Node3DOF}, body,
                           ext_f, du, dr_ids, ω, i, dt, u_t, p, t)
    (a, s) = rod_acceleration(u_v, system, i)
    a = f_acceleration(a, ext_f, i)
    a = constrain_acceleration(a, body)
    a = apply_jns!(a, s, dt)
    return (a, a)
end

function get_system_forces(u_v, system::StructuralGraphSystem{Node6DOF},
                           simulation::RodSimulation{Node6DOF}, body,
                           ext_f, du, dr_ids, ω, i, dt, u_t, p)
    (a, τ, s, j) = rod_acceleration(u_v, system, body, i)

    (a, τ) = f_acceleration(a, τ, ext_f, i, p)
    (a, τ) = constrain_acceleration(a, τ, body)
    s = 0.0 .* s
    j = 0.0 .* j
    a = apply_jns!(a, s, dt)
    dω = update_dω(i, ω, τ, u_v, du, dr_ids, j, u_t, dt)
    return a, dω
end

function get_odeproblem_no_s(simulation::S, ext_f) where {T, S <: StructuralSimulation{T}}
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

        # Set translation vels
        @views du[dx_ids] .= u[v_ids]

        # Get angular velocities
        ω = @view u[ω_ids]

        t_mult = 1.0

        @inbounds for i in 1:n

            # Get current body
            body = bodies[i]
            if typeof(p) == typeof(t_mult)
                t_mult *= p
            end

            # Accelerate system
            (a, dω) = get_system_forces(u_v, system, simulation, body, ext_f, du,
                                        dr_ids, ω, i, dt, u_t, t_mult)
            # Update accelerations
            update_accelerations!(du, a, dω, u_len, i, simulation)
        end
    end

    return ODEProblem(ode_system!, uv0, simulation.tspan)
end