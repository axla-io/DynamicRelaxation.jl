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
    @time sol = solve(prob, alg, dt = simulation.dt, maxiters = maxiters, callback = cb,
                      verbose = false)

    # Plot final state
    u_final = get_state(sol.u[end], u_len, simulation)
    # plot(u_final[1, :], u_final[3, :])

    u_final_true = [0.0 0.5881664318361151 1.1763586829010184 1.7645733133816837 2.352806882286542 2.9410559476265092 3.529317066596153 4.117586795754975 4.705861691208798 5.294138308791203 5.882413204245024 6.470682933403846 7.05894405237349 7.647193117713461 8.235426686618316 8.82364131709898 9.411833568163884 10.0;
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                    0.0 -0.01138579657230201 -0.021348805908338132 -0.02988885311393739 -0.03700578825128829 -0.04269948635209497 -0.04696984742855641 -0.04981679648213399 -0.05124028351014346 -0.05124028351014354 -0.049816796482135495 -0.04696984742855911 -0.042699486352098304 -0.03700578825128928 -0.029888853113940614 -0.02134880590834383 -0.011385796572305722 0.0]
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
