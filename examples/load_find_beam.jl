using Flux, Optim, DiffEqFlux, DiffEqSensitivity


# GENERATE DATA
# -----------------------------------------------------
include("../src/include_lib.jl")

# Define a simple graph system
n_elem = 17
n_pt = n_elem + 1
graph = StaticGraph(path_graph(n_pt))
system = default_system(graph, Node6DOF, :catenary)


# Set loads
#ext_f = point_loads([Pz(-10, system) ], [n_pt], system)
ext_f = uniform_load(Pz(-10_000, system), system)

# Set parameters
maxiters = 500
dt = 0.01
tspan = (0.0, 10.0)

# Create problem
p_true = [1.0]
simulation = LoadScaleRodSimulation{StructuralGraphSystem{Node6DOF},Float64,eltype(ext_f)}(system, tspan, dt, ext_f)
prob = ODEProblem(simulation, p_true)

# Create callback TODO: find a better way
c = 0.7
(_u0, _v0, n, u_len, v_len) = gather_bodies_initial_coordinates(simulation)
(dx_ids, dr_ids, v_ids, ω_ids) = get_vel_ids(u_len, v_len)
affect!(integrator) = affect!(integrator, v_ids, c)
cb = PeriodicCallback(affect!, 1 * dt; initial_affect=true)

# Set algorithm for solver
#alg = Rosenbrock23(autodiff=true)
alg = RK4()

# Solve problem, corresponding to p = 1.0
@time sol = solve(prob, alg, dt=simulation.dt, maxiters=maxiters, callback=cb);
#@profview solve(prob, alg, dt = simulation.dt, maxiters=maxiters, callback = cb);

# Extract final state
u_final = get_state(sol.u[end], u_len)

# Extract initial state (maybe a bit hacky)
u0 = sol.u[1]

# SET UP OPTIMIZATION
# -----------------------------------------------------

# Create a solution (prediction) for a given starting point u0 and set of
# parameters p
function predict(p)
    return get_state(concrete_solve(prob, alg, u0, p, dt=simulation.dt, maxiters=maxiters, callback=cb).u[end], u_len)
end

# Create loss function
function l2loss(p)
    prediction = predict(p)
    loss = sum(abs, prediction .- u_final)
    return loss, prediction
end

function l1loss(p)
    prediction = predict(p)
    loss = sum(abs, prediction .- u_final)
    return loss, prediction
end

# Initial guess of p
p = [0.3]
sol_pred_init = solve(remake(prob, p=p), alg, dt=simulation.dt, maxiters=maxiters, callback=cb)
# Plot final state
u_pred_init = get_state(sol_pred_init.u[end], u_len)


# Callback function to observe training
list_plots = []
iter = 0

callback = function (p, l, pred)
    global iter
    global list_plots

    if iter == 0
        list_plots = []
    end
    iter += 1

    display(l)

    # using `remake` to re-create our `prob` with current parameters `p`
    sol_pred = solve(remake(prob, p=p), alg, dt=simulation.dt, maxiters=maxiters, callback=cb)
    # Plot final state
    u_pred = get_state(sol_pred.u[end], u_len)
    plt = plot(u_final[1, :], u_final[3, :], lw = 1.5, label="Ground Truth")
    plot!(plt, u_pred_init[1, :], u_pred_init[3, :], lw = 1.5, label="Initial Prediction")
    plot!(plt, u_pred[1, :], u_pred[3, :], lw = 1.5, label="Prediction, iter. $iter")
    p_pred = round(p[1], digits = 4)
    plot!(plt, zlims = (-0.3, 0.0), title = "\nLoad Finding, \$p_{\\rm{true}}\$ = 1.0, \$p_{\\rm{pred}}\$ = $p_pred")

    push!(list_plots, plt)


    # Tell sciml_train to not halt the optimization. If return true, then
    # optimization stops.
    return false
end


# OPTIMIZE!
# -----------------------------------------------------

#= result_ode = DiffEqFlux.sciml_train(loss, p,
    BFGS(initial_stepnorm=0.0001), maxiters=10) =#

result_ode = DiffEqFlux.sciml_train(l2loss, p,
    Adam(0.03), maxiters=20, cb=callback)

result_ode = DiffEqFlux.sciml_train(l1loss, result_ode.minimizer,
    BFGS(initial_stepnorm=0.0001), maxiters=10, cb=callback)


# VISUALIZE
# -----------------------------------------------------
p_1 = result_ode.minimizer
loss(p_1)
# Solve problem
@time sol_pred = solve(remake(prob, p=p_1), alg, dt=simulation.dt, maxiters=maxiters, callback=cb);

# Plot final state
u_pred = get_state(sol_pred.u[end], u_len)
plot(u_final[1, :], u_final[3, :], label="Ground Truth")
plot!(u_pred[1, :], u_pred[3, :], label="Prediction")

animate(list_plots, "load_finding.gif",fps=2)
