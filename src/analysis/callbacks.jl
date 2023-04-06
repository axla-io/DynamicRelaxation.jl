
function velocitydecay!(integrator, n, c)
    @views integrator.u[:, (n+1):(2n)] .*= c
end

function velocitydecay!(integrator, n::AbstractVector, c)
    @views integrator.u[n] .*= c
end

function velocityreset!(integrator, n::AbstractVector, c)
    @views integrator.u[n] .*= 0.0
end

# TODO: Make this work
function ke_condition(u, t, integrator, tol, n::AbstractVector)
    c = sum(abs2, (u[n] .- integrator.uprev[n]))
    return c
end