struct BigonRodSimulation <: StructuralSimulation{Node6DOF}
    system::StructuralGraphSystem{Node6DOF}
    tspan::Tuple{Float64, Float64}
    dt::Float64
    end_pts::Vector{Int64}
    k::Float64
    nu::Float64
end

function bigon_forces(a, u, i1, i2, du, simulation)
    x1 = @view u[(7 * (i1 - 1) + 1):(7 * (i1 - 1) + 3)]
    x2 = @view u[(7 * (i2 - 1) + 1):(7 * (i2 - 1) + 3)]

    v1 = @view du[(7 * (i1 - 1) + 1):(7 * (i1 - 1) + 3)]
    v2 = @view du[(7 * (i2 - 1) + 1):(7 * (i2 - 1) + 3)]

    # Calculate contact forces and add to acceleration
    end_distance_vector = x2 - x1
    elastic_force = simulation.k * end_distance_vector
    relative_velocity = v2 - v1
    damping_force = simulation.nu * relative_velocity
    contact_force = elastic_force + damping_force
    a = a + contact_force

    return a
end

function accelerate_system(u_v, system::StructuralGraphSystem{Node6DOF},
                           simulation::BigonRodSimulation, body,
                           ext_f, du, dr_ids, ω, i, dt, u_t, p)
    (a, τ, s, j) = rod_acceleration(u_v, system, body, i)

    (a, τ) = f_acceleration(a, τ, ext_f, i, p)
    (a, τ) = constrain_acceleration(a, τ, body)

    # Apply bigon forces
    if i == simulation.end_pts[1]
        i2 = simulation.end_pts[2]
        u_len = 7 * length(system.bodies)
        a = bigon_forces(a, u_v, i, i2, du, simulation)
    elseif i == simulation.end_pts[2]
        i2 = simulation.end_pts[1]
        u_len = 7 * length(system.bodies)
        a = bigon_forces(a, u_v, i, i2, du, simulation)
    end

    # TODO: Apply bigon torques

    a = apply_jns!(a, s, dt)
    dω = update_dω(i, ω, τ, u_v, du, dr_ids, j, u_t, dt)
    return a, dω
end