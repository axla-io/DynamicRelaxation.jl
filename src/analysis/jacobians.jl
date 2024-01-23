function get_ode_jac(f, u_len, uv0, simulation, jac_p)
	sd = JacPrototypeSparsityDetection(; jac_prototype = jac_p)
	adtype = AutoSparseForwardDiff()
    y = zero(uv0)
    cache = sparse_jacobian_cache(adtype, sd, (dx, x) -> f(dx, x, p, t), y, uv0)
	return (J, u, p, t) -> inner_jac!(J, u, p, t, f, constrained_dofs(u_len, simulation), adtype, cache, y)
end

function inner_jac!(J, u, p, t, f, constrained, adtype,cache, y)
	T = eltype(J)
    sparse_jacobian!(J, adtype, cache, (dx, x) -> f(dx, x, p, t), y, u)
	for k in constrained
		J[k, :] .= T(0.0)
		J[k, k] = T(1.0)
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

function get_jac_prototype(system, u_len, v_len, simulation, apply_constraints)
	dim = u_len + v_len
	jac_p = zeros(dim, dim)
	n_pt = length(system.bodies)

	# Graph independent parts
	for i ∈ 1:n_pt
		x_0 = (i - 1) * 7
		v_0 = u_len + (i - 1) * 6
		# Set positional
		for j in 1:3
			jac_p[x_0+j, v_0+j] = 1.0
		end

		# Set rotational
		jac_p[x_0+4:x_0+7, x_0+1:x_0+4] .= 1.0
		jac_p[x_0+4:x_0+7, v_0+4:v_0+6] .= 1.0
		jac_p[v_0+4:v_0+6, v_0+4:v_0+6] .= 1.0

	end

	# Contribution to dv
	for edge in edges(system.graph)
		a = Int(edge.src)
		b = Int(edge.dst)

		xa_0 = (a - 1) * 7
		xb_0 = (b - 1) * 7
		va_0 = u_len + (a - 1) * 6
		vb_0 = u_len + (b - 1) * 6

		# Contribution to dv
		# b to a
		jac_p[va_0+1:va_0+6, xa_0+1:xa_0+7] .= 1.0
		jac_p[va_0+1:va_0+6, xb_0+1:xb_0+7] .= 1.0

		# a to b
		jac_p[vb_0+1:vb_0+6, xb_0+1:xb_0+7] .= 1.0
		jac_p[vb_0+1:vb_0+6, xa_0+1:xa_0+7] .= 1.0
	end

	# Add constraints
	if apply_constraints
		constrained = constrained_dofs(u_len, simulation)
		for k in constrained
			jac_p[k, :] .= 0.0
			jac_p[k, k] = 1.0
		end
	end

	return sparse(jac_p)
end
