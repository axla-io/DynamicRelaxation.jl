using DynamicRelaxation

using Graphs
using StaticGraphs
using DiffEqCallbacks
using NBodySimulator

using Plots, GraphRecipes
using Test

@testset "Create graph system" begin

    # Define a simple graph and plot
    n_elem = 17
    n_pt = n_elem + 1
    graph = StaticGraph(path_graph(n_pt))

    plot(path_graph(n_pt), curves = false)

    @test nv(graph) == n_pt
    @test ne(graph) == n_elem

end
