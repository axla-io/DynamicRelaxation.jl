using Plots, GraphRecipes

using DynamicRelaxation

using Graphs
using StaticGraphs
using DiffEqCallbacks
using NBodySimulator
using SteadyStateDiffEq
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
simulation = RodSimulation(system, tspan, dt)
prob = ODEProblem(simulation, ext_f)

ssprob = SteadyStateProblem(prob)

# Create callback for damping
c = 0.7
(u0, v0, n, u_len, v_len) = get_u0(simulation)
(dx_ids, dr_ids, v_ids, ω_ids) = get_vel_ids(u_len, v_len, system)
v_decay!(integrator) = velocitydecay!(integrator, vcat(v_ids, ω_ids), c)
cb1 = PeriodicCallback(v_decay!, 2 * dt; initial_affect = true)

# Create callback for termination
abstol = reltol = 1e-2
ktc(integrator, abstol, reltol, min_t) = ke_termination_cond(integrator, abstol, reltol, min_t, v_ids)
cb2 = TerminateSteadyState(abstol, reltol, ktc, min_t = 20*dt)

cb = CallbackSet(cb1,cb2)

# Set algorithm for solver
alg = RK4()

# Solve problem
@time sol = solve(ssprob, DynamicSS(alg), custom_termination_cond = true, maxiters=maxiters, callback = cb);

# Plot final state
u_final = get_state(sol.u, u_len, simulation)
plot(u_final[1, :], u_final[3, :])
