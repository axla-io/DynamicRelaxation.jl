
function velocitydecay!(integrator, n, c)
    @views integrator.u[:, (n+1):(2n)] .*= c
    #integrator.u[:, :] = integrator.u[:, :] * 0.10
end

function velocitydecay!(integrator, n::AbstractVector, c)
    @views integrator.u[n] .*= c
    #integrator.u[:, :] = integrator.u[:, :] * 0.10
end

function velocityreset!(integrator, n::AbstractVector, c)
    @show t
    @views integrator.u[n] .*= 0.0
    #integrator.u[:, :] = integrator.u[:, :] * 0.10
end


function ke_condition(u, t, integrator, tol, n::AbstractVector)
    c = sum(abs2, (u[n] .- integrator.uprev[n])) #./ (t - integrator.tprev) .- tol)
    #@show c
    return c
end

