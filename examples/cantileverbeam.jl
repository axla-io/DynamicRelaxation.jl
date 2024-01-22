using DynamicRelaxation
using OrdinaryDiffEq
using Graphs
using StaticGraphs
using DiffEqCallbacks
using Plots
using LinearAlgebra

n_elem = 6
n_pt = n_elem + 1
graph = StaticGraph(path_graph(n_pt))
system = default_system(graph, Node6DOF, :cantilever, n_pt)

# Set loads
ext_f = point_loads([Pz(-1_000, system)], [n_pt], system)
# Set parameters
maxiters = 5000
dt = 10
tspan = (0.0, Inf)

# Create problem
simulation = RodSimulation(system, tspan, dt)
prob = ODEProblem(simulation, ext_f)

# Create callback
c = 0.991
(u0, v0, n, u_len, v_len) = get_u0(simulation)
(dx_ids, dr_ids, v_ids, ω_ids) = get_vel_ids(u_len, v_len, system)
v_decay!(integrator) = velocitydecay!(integrator, vcat(v_ids, ω_ids), c)
cb1 = PeriodicCallback(v_decay!, 2 * dt; initial_affect = true)

res = zero(prob.u0)
prob.f(res, prob.u0, prob.p, 1.0)
norm(res)
# Set algorithm for solver
alg = RK4()

# Solve problem
@time sol = solve(prob, alg, dt = simulation.dt, maxiters = maxiters, callback = cb1, 
                  verbose = false);
sol.retcode
# Plot final state
u_final = get_state(sol.u[end], u_len, simulation)

# Select frames for animation
itt = generate_range(100, 1, length(sol.u))
u_red = sol.u[itt]

gr()
# Loop over the time values and create a plot for each frame
anim = @animate for i in axes(u_red, 1)
    u_final = get_state(u_red[i], u_len, simulation)
    plot(u_final[1, :], u_final[3, :], xlims = (0,10), ylims = (-5, 1), label = "",c = :black, lw = 3.0)
    quiver!([u_final[1, end]], [u_final[3, end]], quiver=([0], [-1]), color=:blue, label="Arrow Tip", lw = 2.5)
    scatter!(u_final[1, 2:end], u_final[3, 2:end], label = "",markercolor=:white, markerstrokecolor=:black, markersize=5, lw = 3.0)
end

res = zero(prob.u0)
prob.f(res, sol.u[end], prob.p, 0.0)
norm(res)

gif(anim, "cantilever.gif")

