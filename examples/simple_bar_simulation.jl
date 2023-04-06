using DynamicRelaxation

using Graphs
using StaticGraphs
using DiffEqCallbacks
using DifferentialEquations

using Plots, GraphRecipes

# Define a simple graph system
n_elem = 17
n_pt = n_elem + 1
graph = StaticGraph(path_graph(n_pt))
system = default_system(graph, Node3DOF, :catenary, n_pt)

# Set loads
ext_f = uniform_load(Pz(-10, system), system)

# Set parameters
maxiters = 500
dt = 0.001
tspan = (0.0, 10.0)

# Create callback
c = 0.9
v_decay!(integrator) = velocitydecay!(integrator, n_pt, c)
cb = PeriodicCallback(v_decay!,  dt; initial_affect = true)

# Set algorithm for solver
alg = RK4()

# Create problem
simulation = RodSimulation{StructuralGraphSystem{Node3DOF},Float64,eltype(ext_f)}(system, tspan, dt, ext_f)
prob = ODEProblem(simulation)

# Solve problem
@time sol = solve(prob, alg, dt = simulation.dt, maxiters=maxiters, callback = cb);

# Plot final state
u_final = sol.u[end][:, 1:n_pt]
plot(u_final[1, :], u_final[3, :])

# Select frames for animation
itt = generate_range(100, 1, length(sol.u))
u_red = sol.u[itt]

# Loop over the time values and create a plot for each frame
anim = @animate for i in axes(u_red, 1)
    u_final = u_red[i][:, 1:n_pt]
    plot(u_final[1, :], u_final[3, :])
end

# Save the frames as a gif
gif(anim, "animation.gif", fps=20)
