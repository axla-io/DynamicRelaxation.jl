
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
end 

function ke_condition(u, t, integrator, tol, n::AbstractVector)
    c = sum(abs2, u[n]) < sum(abs2, integrator.uprev[n]) && sum(abs2, u[n])>tol
    return c
end
=#

function ke_termination_cond(integrator, abstol, reltol, min_t, n)
    if min_t === nothing
        return sum(abs2, testval)/2 < abstol
    else
        return integrator.iter >= min_t && sum(abs2, integrator.u[n])/2 < abstol
    end
end 

#= function ke_termination_cond(integrator, abstol, reltol, min_t, n)
    if min_t !== nothing
        return integrator.t >= min_t && sum(abs2, get_du(integrator)) < abstol
    end
end
 =#
