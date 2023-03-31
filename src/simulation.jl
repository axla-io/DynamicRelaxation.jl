
abstract type StructureSimulation end

struct RodSimulation{sType<:AbstractGraphSystem,tType<:Real,fType<:Real} <: StructureSimulation
    system::sType
    tspan::Tuple{tType,tType}
    ext_f::Vector{SVector{3,fType}}
end


function gather_bodies_initial_coordinates(simulation::RodSimulation)
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


function DiffEqBase.ODEProblem(simulation::RodSimulation{<:AbstractGraphSystem})
    (u0, v0, n) = gather_bodies_initial_coordinates(simulation)
    bodies = simulation.system.bodies
    system = simulation.system
    ext_f = simulation.ext_f

    function ode_system!(du, u, p, t)
        du[:, 1:n] = @view u[:, (n+1):(2n)]
        #du[:, 1:n] = u[:, (n+1):(2n)]
        #u[:, (n+1):(2n)] .*= 0.99
        u_v = @view u[:, 1:n]
        @inbounds for i in 1:n
            a = MVector(0.0, 0.0, 0.0)
            body = bodies[i]
            #@infiltrate
            rod_acceleration!(a, u_v, system, i)
            f_acceleration!(a, ext_f, i)
            constrain_acceleration!(a, body)
            a = body.M_inv * a

            du[:, n+i] .= a
        end
    end

    return ODEProblem(ode_system!, hcat(u0, v0), simulation.tspan)
end