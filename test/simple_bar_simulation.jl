using DynamicRelaxation

using Graphs
using StaticGraphs
using DiffEqCallbacks
using NBodySimulator

#using Plots, GraphRecipes
using Test

@testset "Bar simulation" begin

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
    cb = PeriodicCallback(v_decay!, dt; initial_affect = true)

    # Set algorithm for solver 
    alg = RK4()

    # Create problem
    simulation = RodSimulation(system, tspan, dt)
    (u0, v0, n, u_len, v_len) = get_u0(simulation)
    prob = ODEProblem(simulation, ext_f)

    # Solve problem
    @time sol = solve(prob, alg, dt = simulation.dt, maxiters = maxiters, callback = cb, verbose = false)

    # Plot final state
    u_final = get_state(sol.u[end], u_len, simulation)
    #plot(u_final[1, :], u_final[3, :])

    u_final_true = [0.0 0.9999178077243016 1.99986643409378 2.9998417722609956 3.999839714390765 4.999856151811957 5.999886975169358 6.9999280745756245 7.999975339763273 9.000024660236727 10.000071925424376 11.000113024830641 12.000143848188042 13.000160285609235 14.000158227739007 15.00013356590622 16.000082192275695 17.0;
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                    0.0 -0.016218327579968574 -0.030409801596914338 -0.042574247125090364 -0.052711514208112986 -0.06082147786820732 -0.06690403811390315 -0.07095911994620795 -0.07298667336321966 -0.07298667336321886 -0.07095911994620602 -0.06690403811389987 -0.06082147786820441 -0.052711514208110134 -0.04257424712508335 -0.03040980159690485 -0.016218327579965747 0.0]

    @test isapprox(sum(abs, u_final .- u_final_true), 0.0, atol = 1e-7)
    #= 
        # Select frames for animation
        itt = generate_range(100, 1, (length(sol.u)-1)/7)
        u_red = sol.u[itt]

        # Loop over the time values and create a plot for each frame
        anim = @animate for i in axes(u_red, 1)
            u_final = get_state(u_red[i], u_len, simulation)
            plot(u_final[1, :], u_final[3, :], label = "", ylims = (-0.06, 0.0))
        end

        # Save the frames as a gif
        gif(anim, "hanging_chain.gif", fps = 20) 
        =#
end
