function get_ode_jac(f, u_len, simulation)
	return (J, u, p, t) -> inner_jac!(J, u, p, t, f, constrained_dofs(u_len, simulation))
end

function inner_jac!(J, u, p, t, f, constrained)
	forwarddiff_color_jacobian!(J, (dx, x) -> f(dx, x, p, t), u)
	for k in constrained
		J[k, :] .= 0.0
		J[k, k] = 1.0
	end
	return nothing
end

function constrained_dofs(u_len, simulation::T) where {T <: StructuralSimulation{Node6DOF}}
	num_constrained_pos = sum(sum(@view(body.constraints[1:3])) for body in simulation.system.bodies if body.constrained == true) * 2 # Times 2 because we constrain velocity
	num_constrained_rot = 7 * (sum(sum(@view(body.constraints[4:7])) for body in simulation.system.bodies if body.constrained == true) / 4) |> Int # Will always either be clamped or free, atleast for now.
	num_constrained = num_constrained_pos + num_constrained_rot
    dofs = zeros(Int, num_constrained)
	ctr = 1
	for (i, body) in enumerate(simulation.system.bodies) # Loop through bodies
		if body.constrained == true
			for (j, cons) in enumerate(body.constraints) # Loop through constraints if the body is constrained
				if cons #if constrained
					dofs[ctr] = (i - 1) * 7 + j # Add position dof to constrained dof list
					ctr += 1
					if j <= 6 # If the dof number is smaller than 6, add the velocity dof. This is because we use rotation vectors for the derivatives.
						dofs[ctr] = u_len + (i - 1) * 6 + j
                        ctr += 1
					end
				end
			end
		end
	end
	return dofs
end


function constrained_dofs(constrained_ids, u_len, simulation::T) where {T <: StructuralSimulation{Node3DOF}}

	return error("constrained_dofs for 3DOF systems not implemented yet!")
end
