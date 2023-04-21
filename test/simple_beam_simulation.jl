using DynamicRelaxation

using Graphs
using StaticGraphs
using DiffEqCallbacks
using NBodySimulator

# using Plots, GraphRecipes
using Test

@testset "Beam simulation" begin

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

    # Create callback
    c = 0.7
    (u0, v0, n, u_len, v_len) = get_u0(simulation)
    (dx_ids, dr_ids, v_ids, ω_ids) = get_vel_ids(u_len, v_len, system)
    v_decay!(integrator) = velocitydecay!(integrator, vcat(v_ids, ω_ids), c)
    cb1 = PeriodicCallback(v_decay!, 2 * dt; initial_affect = true)

    # Set algorithm for solver
    alg = RK4()

    # Solve problem
    @time sol = solve(prob, alg, dt = simulation.dt, maxiters = maxiters, callback = cb1, verbose = false);
    #@profview sol = solve(prob, alg, dt = simulation.dt, maxiters = maxiters, callback = cb1, verbose = false);

    # Plot final state
    u_final = get_state(sol.u[end], u_len, simulation)
    #=
    plot(u_final[1, :], u_final[3, :])
    =#

    u_final_true = [0.0 0.9985840905257405 1.9976994903128098 2.997274873545508 3.997239315074668 4.9975225625937725 5.9980536618600535 6.998761291272035 7.999575190903782 9.000424809096222 10.001238708727966 11.001946338139948 12.002477437406226 13.002760684925333 14.002725126454491 15.002300509687192 16.00141590947426 17.0;
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                    0.0 -0.06146492984809863 -0.1152539050247191 -0.16136341273093313 -0.19979091083461925 -0.23053544160821957 -0.25359487860921576 -0.26896819122986326 -0.27665478826820794 -0.27665478826819967 -0.26896819122984394 -0.25359487860923163 -0.23053544160822184 -0.19979091083461728 -0.16136341273093896 -0.11525390502471794 -0.06146492984809655 0.0]

    @test isapprox(sum(abs, u_final .- u_final_true), 0.0, atol = 1e-4)
    #= 
        # Select frames for animation
        itt = generate_range(100, 1, length(sol.u))
        u_red = sol.u[itt]

        # Loop over the time values and create a plot for each frame
        anim = @animate for i in axes(u_red, 1)
            u_final = get_state(u_red[i], u_len, simulation)
            plot(u_final[1, :], u_final[3, :])
        end

        # Save the frames as a gif
        gif(anim, fps = 20) 
        =#
end
