using Optimization, SciMLSensitivity, Zygote

using OptimizationOptimJL, OptimizationOptimisers

using Plots, GraphRecipes

using DynamicRelaxation

using Graphs
using StaticGraphs
using DiffEqCallbacks
using NBodySimulator
using DiffEqFlux

# GENERATE DATA
# -----------------------------------------------------

# Define a simple graph system
n_elem = 17
n_pt = n_elem + 1
graph = StaticGraph(path_graph(n_pt))
system = default_system(graph, Node6DOF, :catenary, n_pt)

# Set loads
ext_f = uniform_load(Pz(-10_000, system), system)

# Set parameters
maxiters = 500
dt = 0.01
tspan = (0.0, 10.0)

# Create problem
simulation = LoadScaleRodSimulation(system, tspan, dt)
prob = ODEProblem(simulation, ext_f)

# Create decay callback
c = 0.7
(_u0, _v0, n, u_len, v_len) = get_u0(simulation)
(dx_ids, dr_ids, v_ids, ω_ids) = get_vel_ids(u_len, v_len, system)
v_decay!(integrator) = velocitydecay!(integrator, v_ids, c)

cb1 = PeriodicCallback(v_decay!, 1 * dt; initial_affect = true)

# Set algorithm for solver
alg = RK4()

# Solve problem, corresponding to p = 1.0
p_true = [1.0]
@time sol = solve(prob, p = p_true, alg, dt = simulation.dt, maxiters = maxiters);#, callback=cb1);
# Extract final state
u_final = get_state(sol.u[end], u_len, simulation)

# Extract initial state
u0 = sol.u[1]

# SET UP OPTIMIZATION
# -----------------------------------------------------

# Create a solution (prediction) for a given starting point u0 and set of
# parameters p
#= function predict(p)
    return solve(prob, alg, p = p, dt = simulation.dt, maxiters = maxiters, callback = cb)
end =#

# Create loss function
function l2loss(p)
    prediction = solve(prob, alg, p = p, dt = simulation.dt, maxiters = maxiters)#, callback = cb1)
    u_pred = get_state(prediction.u[end], u_len, simulation)
    loss = sum(abs2, u_pred .- u_final)
    return loss, prediction
end

function l1loss(p)
    prediction = solve(prob, alg, p = p, dt = simulation.dt, maxiters = maxiters)#, callback = cb1)
    u_pred = get_state(prediction.u[end], u_len, simulation)
    loss = sum(abs, u_pred .- u_final)
    return loss, prediction
end

# Initial guess of p
p0 = [0.3]
sol_pred_init = solve(prob, alg, p = p0, dt = simulation.dt, maxiters = maxiters)#,
#callback = cb1);

l2loss(p0)
# Plot final state
u_pred_init = get_state(sol_pred_init.u[end], u_len, simulation)

# Callback function to observe training
list_plots = []
iter = 0

callback = function (p, l, prediction)
    global iter
    global list_plots

    if iter == 0
        list_plots = []
    end
    iter += 1

    display(l)
    u_pred = get_state(prediction.u[end], u_len, simulation)

    # using `remake` to re-create our `prob` with current parameters `p`
    #sol_pred = solve(remake(prob, p=p), alg, dt=simulation.dt, maxiters=maxiters, callback=cb)
    # Plot final state
    #u_pred = get_state(sol_pred.u[end], u_len)
    plt = plot(u_final[1, :], u_final[3, :], lw = 1.5, label = "Ground Truth")
    plot!(plt, u_pred_init[1, :], u_pred_init[3, :], lw = 1.5, label = "Initial Prediction")
    plot!(plt, u_pred[1, :], u_pred[3, :], lw = 1.5, label = "Prediction, iter. $iter")
    p_pred = round(p[1], digits = 4)
    plot!(plt, ylims = (-0.3, 0.0),
          title = "\nLoad Finding, \$p_{\\rm{true}}\$ = 1.0, \$p_{\\rm{pred}}\$ = $p_pred")

    push!(list_plots, plt)

    # Tell sciml_train to not halt the optimization. If return true, then
    # optimization stops.
    return false
end

# OPTIMIZE!
# -----------------------------------------------------

function axl_train(loss, θ, opt = OptimizationPolyalgorithms.PolyOpt(), adtype = nothing,
                   args...;
                   lower_bounds = nothing, upper_bounds = nothing, cb = nothing,
                   callback = (args...) -> (false),
                   maxiters = nothing, kwargs...)
    adtype = Optimization.AutoFiniteDiff()
    #adtype = Optimization.AutoZygote()

    if !isnothing(cb)
        callback = cb
    end

    optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
    optprob = Optimization.OptimizationProblem(optf, θ; lb = lower_bounds,
                                               ub = upper_bounds, kwargs...)

    Optimization.solve(optprob, opt, args...; maxiters, callback = callback, kwargs...)
end

result_ode1 = axl_train(l2loss, p0,
OptimizationOptimisers.Adam(0.03), maxiters=20)

#= result_ode1 = DiffEqFlux.sciml_train(l2loss, p0,
OptimizationOptimisers.Adam(0.03), maxiters=20) =#
#= 
adtype = Optimization.AutoZygote()
optf1 = Optimization.OptimizationFunction((x, p) -> l2loss(x), adtype)
optprob1 = Optimization.OptimizationProblem(optf1, p0)

result_ode1 = Optimization.solve(optprob1, Adam(0.03),
                                 #callback = callback,
                                 maxiters = 20)

                               @report_call  Optimization.solve(optprob1, Adam(0.03),
                                 #callback = callback,
                                 maxiters = 20) =#
#= #Remake and solve with BFGS
optf2 = Optimization.OptimizationFunction((x, p) -> l1loss(x), adtype)
optprob2 = Optimization.OptimizationProblem(optf2, result_ode1.minimizer)

result_ode2 = Optimization.solve(optprob2, BFGS(initial_stepnorm = 0.0001),
                                 callback = callback,
                                 maxiters = 10)
 =#
# VISUALIZE
# -----------------------------------------------------
p_1 = result_ode1.minimizer

# Solve problem
@time sol_pred = solve(remake(prob, p = p_1), alg, dt = simulation.dt, maxiters = maxiters,
                       callback = cb);

# Plot final state
u_pred = get_state(sol_pred.u[end], u_len, simulation)
plot(u_final[1, :], u_final[3, :], label = "Ground Truth")
plot!(u_pred[1, :], u_pred[3, :], label = "Prediction")

animate(list_plots, "load_finding.gif", fps = 2)
