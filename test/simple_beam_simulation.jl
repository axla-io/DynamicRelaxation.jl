using DynamicRelaxation

using Graphs
using StaticGraphs
using DiffEqCallbacks
using NBodySimulator

#using Plots, GraphRecipes
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
    cb1 = PeriodicCallback(v_decay!, 3 * dt; initial_affect = true)

    # Set algorithm for solver
    alg = RK4()

    # Solve problem
    @time sol = solve(prob, alg, dt = simulation.dt, maxiters = maxiters, callback = cb1)

    # Plot final state
    u_final = get_state(sol.u[end], u_len, simulation)
    #plot(u_final[1, :], u_final[3, :])

    u_final_true = [0.0 0.9986722089024265 1.9977858647654367 2.997273043411182 3.997143682802792 4.997451957090221 5.998094369849439 6.998858956587347 7.999628750605406 9.000371249394592 10.00114104341265 11.001905630150565 12.00254804290978 13.002856317197205 14.00272695658882 15.002214135234565 16.001327791097573 17.0;
                    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
                    0.0 -0.061571676604151505 -0.11539880745740831 -0.1612191993541145 -0.199801504907901 -0.23036953682824707 -0.25354819106840504 -0.26928831873876025 -0.27716124078038057 -0.2771612407803893 -0.2692883187387831 -0.2535481910684114 -0.2303695368282427 -0.1998015049078976 -0.16121919935412593 -0.11539880745739693 -0.06157167660414585 0.0]

    @test isapprox(sum(abs, u_final .- u_final_true), 0.0, atol = 1e-7)
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
    gif(anim, fps = 20) =#
end
