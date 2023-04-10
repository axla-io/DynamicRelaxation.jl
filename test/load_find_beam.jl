using Optimization, SciMLSensitivity, Zygote

using OptimizationNLopt

using DynamicRelaxation

using Graphs
using StaticGraphs
using DiffEqCallbacks
using NBodySimulator

#using Plots, GraphRecipes
using Test

# For plotting
list_plots = []
iter = 0

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
maxiters = 1000
dt = 0.01
tspan = (0.0, 10.0)

# Create problem
simulation = LoadScaleRodSimulation(system, tspan, dt)
prob = ODEProblem(simulation, ext_f)

# Create decay callback
c = 0.7
(_u0, _v0, n, u_len, v_len) = get_u0(simulation)
(dx_ids, dr_ids, v_ids, ω_ids) = get_vel_ids(u_len, v_len, system)
v_decay!(integrator) = velocitydecay!(integrator, vcat(v_ids, ω_ids), c)

cb1 = PeriodicCallback(v_decay!, 1 * dt; initial_affect = true)

# Set algorithm for solver
alg = RK4()

# Solve problem, corresponding to p = 1.0
p_true = [1.0]
sol = solve(prob, p = p_true, alg, dt = simulation.dt, maxiters = maxiters,
                  callback = cb1, verbose = false)
# Extract final state
u_final = get_state(sol.u[end], u_len, simulation)

# SET UP OPTIMIZATION
# -----------------------------------------------------

# Create a solution (prediction) for a given starting point u0 and set of
# parameters p
function predict(p)
    return get_state(solve(prob, alg, p = p, dt = simulation.dt, maxiters = maxiters, callback = cb1, verbose = false).u[end],
                     u_len, simulation)
end

# Create loss function
function l2loss(p)
    u_pred = predict(p)
    loss = sum(abs2, u_pred .- u_final)
    return loss, u_pred
end

# Initial guess of p
p0 = [0.3]
sol_pred_init = solve(prob, alg, p = p0, dt = simulation.dt, maxiters = maxiters,
                      callback = cb1, verbose = false)

# Plot final state
u_pred_init = get_state(sol_pred_init.u[end], u_len, simulation)

# Callback function to observe training

#= 
callback = function (p, l, prediction)
    global iter
    global list_plots

    if iter == 0
        list_plots = []
    end
    iter += 1

    display(l)
    
    u_pred = prediction
    plt = plot(u_final[1, :], u_final[3, :], lw = 1.5, label = "Ground Truth")
    plot!(plt, u_pred_init[1, :], u_pred_init[3, :], lw = 1.5,
        label = "Initial Prediction")
    plot!(plt, u_pred[1, :], u_pred[3, :], lw = 1.5, label = "Prediction, iter. $iter")
    p_pred = round(p[1], digits = 4)
    plot!(plt, ylims = (-0.4, 0.0),
        title = "\nLoad Finding, \$p_{\\rm{true}}\$ = 1.0, \$p_{\\rm{pred}}\$ = $p_pred")

    push!(list_plots, plt) 
    

    return false
end
=#

# OPTIMIZE!
# -----------------------------------------------------
adtype = Optimization.AutoForwardDiff()
optf1 = Optimization.OptimizationFunction((x, p) -> l2loss(x), adtype)
optprob1 = Optimization.OptimizationProblem(optf1, p0)

@time result_ode = Optimization.solve(optprob1, NLopt.LD_LBFGS(),
                                      #callback = callback,
                                      maxiters = 6)

# VISUALIZE
# -----------------------------------------------------
p_1 = result_ode.u

# Solve problem
sol_pred = solve(remake(prob, p = p_1), alg, dt = simulation.dt,
                       maxiters = maxiters,
                       callback = cb1, verbose = false)

# Plot final state
u_pred = get_state(sol_pred.u[end], u_len, simulation)

#= 
plot(u_final[1, :], u_final[3, :], label = "Ground Truth")
plot!(u_pred[1, :], u_pred[3, :], label = "Prediction") 
=#

@testset "Load find beam" begin
    @test isapprox(p_1[1], 1.0, atol = 1e-5)
    @test isapprox(l2loss(p_1)[1], 0.0, atol = 1e-5)
    @test isapprox(sum(abs, u_final .- u_pred), 0.0, atol = 1e-4)
end

#= 
animate(list_plots, "load_finding_bfgs.gif", fps = 2) 
=#
