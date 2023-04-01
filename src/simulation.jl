
abstract type StructureSimulation end

struct RodSimulation{sType<:AbstractGraphSystem,tType<:Real,fType<:Real} <: StructureSimulation
    system::sType
    tspan::Tuple{tType,tType}
    dt::tType
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

function generate_range(n, t1, t2)
    step = (t2 - t1) / (n - 1)
    return [round(Int, t1 + i * step) for i in 0:n-1]
end

function condition(u, t, integrator)
    #return (integrator.t>0.3)
    return true
end

function affect!(integrator, n, c)
    @views integrator.u[:, (n+1):(2n)] .*= c
    #integrator.u[:, :] = integrator.u[:, :] * 0.10
end