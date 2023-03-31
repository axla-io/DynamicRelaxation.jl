using Infiltrator

include("simple_graph.jl")


function generate_range(n, t1, t2)
    step = (t2 - t1) / (n - 1)
    return [round(Int, t1 + i*step) for i in 0:n-1]
end

tspan = (0.0, 100.0)
ext_f = [SVector{3,Float64}([0.0, 0.0, -10.0]) for i in 1:Int(nv(system.graph))]

simulation = RodSimulation{StructuralGraphSystem{Node3DOF},Float64,Float64}(system, tspan, ext_f)

prob = ODEProblem(simulation)
system.bodies
#alg = Rosenbrock23(autodiff=false);
#alg = RK4()
#maxiters = 3000000
alg = ImplicitEuler(autodiff=false)
maxiters = 1000
abstol = 1e-1
dt = 0.01

@time sol = solve(prob, alg, dt = dt, maxiters=maxiters);
#@time sol = solve(prob, alg, maxiters=maxiters);
u_final = sol.u[end][:, 1:n_pt]

plot(u_final[1, :], u_final[3, :])
itt = generate_range(100, 1, length(sol.u))
u_red = sol.u[itt]

# Loop over the time values and create a plot for each frame
anim = @animate for i in axes(u_red, 1)
    u_final = u_red[i][:, 1:n_pt]
    plot(u_final[1, :], u_final[3, :])
end

# Save the frames as a gif
gif(anim, "animation.gif", fps=20)