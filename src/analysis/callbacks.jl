
function velocitydecay!(integrator, n, c)
    @views integrator.u[(3n + 1):(6n)] .*= c
end

function velocitydecay!(integrator, n::AbstractVector, c)
    @views integrator.u[n] .*= c
end

function velocityreset!(integrator, n::AbstractVector, c)
    @views integrator.u[n] .*= 0.0
    println("hello")
end

# TODO: Make this work
#= function ke_condition(u, t, integrator, tol, n::AbstractVector)
    c = sum(abs2, u[n]) - sum(abs2, integrator.uprev[n]) - tol
    return c
end =#

function ke_condition(u, t, integrator, tol, n::AbstractVector)
    c = sum(abs2, u[n]) < sum(abs2, integrator.uprev[n]) && sum(abs2, u[n])>tol
    return c
end