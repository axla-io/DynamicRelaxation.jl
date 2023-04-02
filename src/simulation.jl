
abstract type StructureSimulation end

struct RodSimulation{sType<:AbstractGraphSystem,tType<:Real,fType} <: StructureSimulation
    system::sType
    tspan::Tuple{tType,tType}
    dt::tType
    ext_f::Vector{fType}
end

function gather_bodies_initial_coordinates(simulation::RodSimulation{StructuralGraphSystem{Node3DOF},Float64,SVector{3,Float64}})
    system = simulation.system
    bodies = system.bodies
    len = n = length(bodies)

    u0 = zeros(3, len)
    v0 = zeros(3, len)

    for i in 1:n
        u0[:, i] = bodies[i].r
        v0[:, i] = bodies[i].v
    end

    (u0, v0, n)
end


function gather_bodies_initial_coordinates(simulation::RodSimulation{StructuralGraphSystem{Node6DOF},Float64,SVector{6,Float64}})
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
        u0[bx_id:bx_id+2] = bodies[i].r
        u0[bx_id+3:bx_id+6] = bodies[i].q
        v0[bv_id:bv_id+2] = bodies[i].v
        v0[bv_id+3:bv_id+5] = bodies[i].ω
    end

    (u0, v0, n, u_len, v_len)
end


function DiffEqBase.ODEProblem(simulation::RodSimulation{StructuralGraphSystem{Node3DOF},Float64,SVector{3,Float64}})
    (u0, v0, n) = gather_bodies_initial_coordinates(simulation)
    bodies = simulation.system.bodies
    system = simulation.system
    ext_f = simulation.ext_f
    dt = simulation.dt

    function ode_system!(du, u, p, t)
        du[:, 1:n] = @view u[:, (n+1):(2n)]
        u_v = @view u[:, 1:n]
        u_t = eltype(u)
        a = @MVector zeros(u_t, 3)
        s = @MVector zeros(u_t, 3)
        _z = zero(u_t)
        @inbounds for i in 1:n
            @views a .*= _z
            @views s .*= _z
            body = bodies[i]
            #@infiltrate
            rod_acceleration!(a, u_v, system, i, s)
            f_acceleration!(a, ext_f, i)
            constrain_acceleration!(a, body)
            s .*= dt^2.0 / 2.0
            s_min!(s)
            a = a ./ s

            du[:, n+i] .= a
        end
    end

    return ODEProblem(ode_system!, hcat(u0, v0), simulation.tspan)
end

function DiffEqBase.ODEProblem(simulation::RodSimulation{StructuralGraphSystem{Node6DOF},Float64,SVector{6,Float64}})
    (u0, v0, n, u_len, v_len) = gather_bodies_initial_coordinates(simulation)
    uv0 = vcat(u0, v0)
    (dx_ids, dr_ids, v_ids, ω_ids) = get_vel_ids(u_len, v_len)

    bodies = simulation.system.bodies
    system = simulation.system
    ext_f = simulation.ext_f
    dt = simulation.dt

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
        set_rotation_vels!(dr, ω, n)

        # Initialize velocity and acceleration
        a = @MVector zeros(u_t, 3)
        s = @MVector zeros(u_t, 3)
        j = @MVector ones(u_t, 3) # change and also add j update
        τ = @MVector zeros(u_t, 3) # change and also add τ update

        @inbounds for i in 1:n
            # Reset accelerations
            @views a .*= _z
            @views s .*= _z
            dω_id = 3 * (i - 1) + 1
            ω_i = SVector{3,u_t}(ω[dω_id:dω_id+2])

            body = bodies[i]
            rod_acceleration!(a, τ, u_v, system, body, i, s, j)
            f_acceleration!(a, τ, ext_f, i)
            constrain_acceleration!(a, τ, body)

            # Apply momentum and masses
            apply_jns!(a, s, dt)
            apply_jns!(τ, j, dt)

            #@infiltrate
            # Apply moment of inertia
            dω = (τ - j .* scross(ω_i, ω_i)) ./ j

            # Update accelerations
            d_id = (u_len) + 6 * (i - 1) + 1
            @views du[d_id:d_id+2] .= a
            @views du[d_id+3:d_id+5] .= dω

        end
    end

    return ODEProblem(ode_system!, uv0, simulation.tspan)
end


function generate_range(n, t1, t2)
    step = (t2 - t1) / (n - 1)
    return [round(Int, t1 + i * step) for i in 0:n-1]
end

function apply_jns!(a, s, dt)
    s .*= dt^2.0 / 2.0
    s_min!(s)
    a .= a ./ s
    return nothing
end

function condition(u, t, integrator)
    #return (integrator.t>0.3)
    return true
end

function affect!(integrator, n, c)
    @views integrator.u[:, (n+1):(2n)] .*= c
    #integrator.u[:, :] = integrator.u[:, :] * 0.10
end

function affect!(integrator, n::AbstractVector, c)
    @views integrator.u[n] .*= c
    #integrator.u[:, :] = integrator.u[:, :] * 0.10
end

function get_vel_ids(u_len, v_len)

    dx_ids = get_ids(1, 3, 7, u_len)
    dr_ids = get_ids(4, 4, 7, u_len)
    v_ids = get_ids(u_len + 1, 3, 6, u_len + v_len)
    ω_ids = get_ids(u_len + 3, 3, 6, u_len + v_len)

    return dx_ids, dr_ids, v_ids, ω_ids
end

function get_state(u, u_len)
    x_ids = 1:7:u_len
    y_ids = 2:7:u_len
    z_ids = 3:7:u_len

    state = hcat(u[x_ids], u[y_ids], u[z_ids])'

    return state
end

function get_ids(start, step_inc, offset, finish)
    rangelist = Vector{Vector{Int64}}()

    for i in 1:step_inc
        push!(rangelist, collect(start+i-1:offset:finish))
    end
    return hcat(rangelist...)'[:]
end

