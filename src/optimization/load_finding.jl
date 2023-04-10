struct LoadScaleRodSimulation{T}<: StructuralSimulation{T}
    system::StructuralGraphSystem{T}
    tspan::Tuple{Float64, Float64}
    dt::Float64
end

function LoadScaleRodSimulation(system::T, tspan, dt) where {T}
    return LoadScaleRodSimulation{eltype(system.bodies)}(system::T, tspan, dt)
end

function f_acceleration!(a, τ, ext_f, i, p)
    for j = 1:3
        a[j] += p[1] * ext_f[i][j]
        τ[j] += p[1] * ext_f[i][j+3]
    end
    return nothing
end

function accelerate_system!(a, τ, dω, u_v, system::StructuralGraphSystem{Node6DOF}, simulation::LoadScaleRodSimulation{Node6DOF}, body,
    ext_f, dr, ω, i, s, j, dt, u_t, p)
rod_acceleration!(a, τ, u_v, system, body, i, s, j)
f_acceleration!(a, τ, ext_f, i, p)
constrain_acceleration!(a, τ, body)
apply_jns!(a, s, dt)
update_dω!(i, ω, τ, dr, j, dω, u_t, dt)
return nothing
end
