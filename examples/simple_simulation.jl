
using Plots, GraphRecipes

include("../src/include_lib.jl")

# Define a simple graph system
n_elem = 17
n_pt = n_elem + 1
graph = StaticGraph(path_graph(n_pt))
system = default_system(graph)

# Set loads
ext_f = uniform_load([0.0, 0.0, -10.0], system)

# Set parameters
maxiters = 1000
dt = 0.0001
tspan = (0.0, 10.0)

# Create callback
c = 0.995
affect!(integrator) = affect!(integrator, n_pt, c)
cb = PeriodicCallback(affect!, dt; initial_affect = true)

# Set algorithm for solver
alg = TRBDF2(autodiff = false)

# Create problem
simulation = RodSimulation{StructuralGraphSystem{Node3DOF},Float64,Float64}(system, tspan, dt, ext_f)
prob = ODEProblem(simulation)

# Solve problem
@time sol = solve(prob, alg, dt = simulation.dt, maxiters=maxiters, callback = cb);
# @time sol = solve(prob, alg, maxiters=maxiters);

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