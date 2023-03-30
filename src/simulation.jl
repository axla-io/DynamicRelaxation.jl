
abstract type StructureSimulation end

struct RodSimulation{sType <: AbstractGraphSystem, tType <: Real, fType<:Real} <: StructureSimulation
system::sType
tspan::Tuple{tType, tType}
ext_f::Vector{SVector{3, fType}}
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


function DiffEqBase.ODEProblem(simulation<:RodSimulation{<:AbstractGraphSystem})
    (u0, v0, n) = gather_bodies_initial_coordinates(simulation)

    acceleration_functions = gather_accelerations(simulation)

    ode_system! = let acceleration_functions = tuple(acceleration_functions...)
        function ode_system!(du, u, p, t)
            du[:, 1:n] = @view u[:, (n + 1):(2n)]

            @inbounds for i in 1:n
                a = MVector(0.0, 0.0, 0.0)
                for acceleration! in acceleration_functions
                    acceleration!(a, u[:, 1:n], u[:, (n + 1):end], t, i)
                end
                du[:, n + i] .= a
            end
        end
    end

    return ODEProblem(ode_system!, hcat(u0, v0), simulation.tspan)
end