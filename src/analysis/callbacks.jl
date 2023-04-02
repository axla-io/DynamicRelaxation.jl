
function velocitydecay!(integrator, n, c)
    @views integrator.u[:, (n+1):(2n)] .*= c
    #integrator.u[:, :] = integrator.u[:, :] * 0.10
end

function velocitydecay!(integrator, n::AbstractVector, c)
    @views integrator.u[n] .*= c
    #integrator.u[:, :] = integrator.u[:, :] * 0.10
end
