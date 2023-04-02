
struct LoadScaleRodSimulation{sType<:AbstractGraphSystem,tType<:Real,fType} <: StructureSimulation
    system::sType
    tspan::Tuple{tType,tType}
    dt::tType
    ext_f::Vector{fType}
end

function f_acceleration!(a, τ, ext_f, i, p)
    for j = 1:3
        a[j] += p[1] * ext_f[i][j]
        τ[j] += p[1] * ext_f[i][j+3]
    end
    return nothing
end

# Change this!
function gather_bodies_initial_coordinates(simulation::LoadScaleRodSimulation{StructuralGraphSystem{Node6DOF},Float64,SVector{6,Float64}})
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



function DiffEqBase.ODEProblem(simulation::LoadScaleRodSimulation{StructuralGraphSystem{Node6DOF},Float64,SVector{6,Float64}}, p)
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
        τ = @MVector zeros(u_t, 3) # change and also add τ update
        j = @MVector zeros(u_t, 3) # change and also add j update

        @inbounds for i in 1:n
            # Reset accelerations
            @views a .*= _z
            @views s .*= _z
            @views τ .*= _z
            @views j .*= _z
            dω_id = 3 * (i - 1) + 1
            ω_i = SVector{3,u_t}(ω[dω_id:dω_id+2])

            #= if i == 2
                @infiltrate
            end =#

            body = bodies[i]
            rod_acceleration!(a, τ, u_v, system, body, i, s, j)
            f_acceleration!(a, τ, ext_f, i, p)
            constrain_acceleration!(a, τ, body)

            # Apply momentum and masses
            apply_jns!(a, s, dt)

            #@infiltrate
            # Apply moment of inertia
            j .*= dt^2.0 / 2.0
            s_min!(j)
            dω = (τ -  scross(ω_i, j .* ω_i)) ./ j

            # Update accelerations
            d_id = (u_len) + 6 * (i - 1) + 1
            @views du[d_id:d_id+2] .= a
            @views du[d_id+3:d_id+5] .= dω

        end
    end

    return ODEProblem(ode_system!, uv0, simulation.tspan, p)
end


