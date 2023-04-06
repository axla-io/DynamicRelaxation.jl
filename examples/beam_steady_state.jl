using Plots, GraphRecipes

using DynamicRelaxation

using Graphs
using StaticGraphs
using DiffEqCallbacks
using DifferentialEquations
using SteadyStateDiffEq

# Define a simple graph system
n_elem = 17
n_pt = n_elem + 1
graph = StaticGraph(path_graph(n_pt))
system = default_system(graph, Node6DOF, :catenary, n_pt)

# Set loads
ext_f = uniform_load(Pz(-10_000, system), system)

# Set parameters
maxiters = 300
dt = 0.01
tspan = (0.0, 10.0)

# Create problem
simulation = RodSimulation{StructuralGraphSystem{Node6DOF},Float64,eltype(ext_f)}(system, tspan, dt, ext_f)
prob = ODEProblem(simulation)
ssprob = SteadyStateProblem(prob)

# Create decay callback
c = 0.7
(u0, v0, n, u_len, v_len) = get_u0(simulation)
(dx_ids, dr_ids, v_ids, ω_ids) = get_vel_ids(u_len, v_len)
v_decay!(integrator) = velocitydecay!(integrator, vcat(v_ids, ω_ids), c)
cb = PeriodicCallback(v_decay!, 3 * dt; initial_affect=true)

# Set termination condition
cond = NLSolveTerminationCondition(NLSolveTerminationMode.AbsSafe; abstol = 1e-1)

# Set algorithm for solver
alg = RK4()

# Solve problem
@time sol = solve(ssprob, DynamicSS(alg, termination_condition = cond),  maxiters=maxiters, callback = cb);

# Plot final state
u_final = get_state(sol.u, u_len)
plot(u_final[1, :], u_final[3, :])
