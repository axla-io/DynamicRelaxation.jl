struct KETerminationCondition{T}<:TerminationCondition
abstol::T
reltol::T
v_inds::Vector{Int}
end

function KETerminationCondition(v_inds; abstol = 1e-8, reltol = 1e-6)
    return KETerminationCondition{typeof(abstol)}(abstol, reltol, v_inds)
end

DiffEqBase.get_termination_mode(cond::KETerminationCondition) = nothing

function DiffEqBase._has_converged(du, u, uprev, cond::KETerminationCondition) # Hacky
    if norm(u[cond.v_inds]) < cond.abstol
        return true
    else
        return false
    end

end

function (cond::KETerminationCondition)(storage::Union{<:AbstractDict, Nothing})

    function test_f(integrator, abstol, reltol, min_t)
        _u = integrator.u
        #@show(norm(_u[cond.v_inds]))
        if norm(_u[cond.v_inds]) < abstol
            return true
        else
            return false
        end
    
    end

    return test_f
end


